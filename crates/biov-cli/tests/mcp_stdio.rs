//! End-to-end tests launch the actual Rust binary and speak newline-delimited
//! JSON-RPC over its pipes. No Python interpreter or in-process mock is used.
use std::{
    collections::BTreeSet,
    fs,
    io::Cursor,
    path::{Path, PathBuf},
    process::Stdio,
    sync::atomic::{AtomicU64, Ordering},
    time::{Duration, SystemTime, UNIX_EPOCH},
};

use base64::{engine::general_purpose::STANDARD, Engine};
use md5::Md5;
use polars::prelude::*;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use tokio::{
    io::{AsyncBufReadExt, AsyncReadExt, AsyncWriteExt, BufReader},
    process::{Child, ChildStdin, ChildStdout, Command},
    time::timeout,
};

const DEADLINE: Duration = Duration::from_secs(20);
const CSV: &str = "gene,score,signal,active\nBRCA1,5,low,true\nTP53,30,high,false\nEGFR,20,mid,true\nALK,30,high,true\nKRAS,,unknown,false\n";
static NEXT_TEMP: AtomicU64 = AtomicU64::new(1);

// CI can run this same contract suite against an independently installed binary.
fn binary() -> PathBuf {
    std::env::var_os("BIOV_TEST_BINARY")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_BIN_EXE_biov")))
}

struct Fixture {
    root: PathBuf,
    data: PathBuf,
    output: PathBuf,
}
impl Fixture {
    fn new() -> Self {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "biov-mcp-test-{}-{nonce}-{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        let data = root.join("data");
        let output = root.join("output");
        fs::create_dir_all(&data).unwrap();
        fs::create_dir_all(&output).unwrap();
        fs::write(data.join("genes.csv"), CSV).unwrap();
        Self { root, data, output }
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

struct Server {
    child: Child,
    stdin: Option<ChildStdin>,
    stdout: BufReader<ChildStdout>,
    sequence: u64,
}
impl Server {
    async fn start(fixture: &Fixture) -> Self {
        Self::start_roots(&fixture.data, &fixture.output).await
    }
    async fn start_roots(data: &Path, output: &Path) -> Self {
        Self::start_with_storage(data, output, None).await
    }
    async fn start_with_storage(data: &Path, output: &Path, store: Option<&Path>) -> Self {
        let mut command = Command::new(binary());
        command
            .arg("mcp-native")
            .arg("--data-root")
            .arg(data)
            .arg("--output-root")
            .arg(output);
        if let Some(store) = store {
            command.arg("--store-root").arg(store);
        }
        let mut child = command
            // An empty PATH also catches accidental Python/helper-process dependencies.
            .env("PATH", "")
            .env_remove("PYTHONPATH")
            .env("POLARS_MAX_THREADS", "2")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .kill_on_drop(true)
            .spawn()
            .unwrap();
        let stdin = child.stdin.take().unwrap();
        let stdout = BufReader::new(child.stdout.take().unwrap());
        let mut server = Self {
            child,
            stdin: Some(stdin),
            stdout,
            sequence: 0,
        };
        let initialization = server
            .request(
                "initialize",
                json!({
                    "protocolVersion": "2025-06-18", "capabilities": {},
                    "clientInfo": {"name": "biov-rust-integration-test", "version": "1"}
                }),
            )
            .await;
        let result = initialization.get("result").expect("initialize result");
        assert_eq!(result["serverInfo"]["name"], "biov-native");
        assert!(result["capabilities"]["tools"].is_object());
        assert_eq!(result["protocolVersion"], "2025-06-18");
        server
            .send(json!({"jsonrpc":"2.0", "method":"notifications/initialized"}))
            .await;
        server
    }

    async fn send(&mut self, message: Value) {
        let bytes = serde_json::to_vec(&message).unwrap();
        let stdin = self.stdin.as_mut().unwrap();
        timeout(DEADLINE, async {
            stdin.write_all(&bytes).await.unwrap();
            stdin.write_all(b"\n").await.unwrap();
            stdin.flush().await.unwrap();
        })
        .await
        .expect("write timed out");
    }

    async fn request(&mut self, method: &str, params: Value) -> Value {
        self.sequence += 1;
        let id = format!("request-{}", self.sequence);
        self.send(json!({"jsonrpc":"2.0", "id": id, "method": method, "params": params}))
            .await;
        timeout(DEADLINE, async {
            loop {
                let mut line = String::new();
                let count = self.stdout.read_line(&mut line).await.unwrap();
                assert!(count > 0, "MCP stdout closed before responding to {method}");
                let message: Value = serde_json::from_str(&line).unwrap_or_else(|error| {
                    panic!("non-JSON output corrupted MCP stdout: {error}: {line}")
                });
                assert_eq!(message["jsonrpc"], "2.0");
                if message.get("id").is_some() {
                    assert_eq!(message["id"], id, "response ID must be preserved");
                    return message;
                }
                assert!(message["method"].is_string(), "invalid notification");
            }
        })
        .await
        .unwrap_or_else(|_| panic!("MCP response timed out for {method}"))
    }

    async fn call(&mut self, name: &str, arguments: Value) -> Value {
        self.request("tools/call", json!({"name": name, "arguments": arguments}))
            .await
    }

    async fn successful(&mut self, name: &str, arguments: Value) -> Value {
        let response = self.call(name, arguments).await;
        assert!(
            response.get("error").is_none(),
            "protocol error: {response}"
        );
        assert!(
            serde_json::to_vec(&response).unwrap().len() <= 96 * 1024,
            "wire response exceeds budget"
        );
        let result = &response["result"];
        assert_eq!(result["isError"], false, "tool failed: {response}");
        assert!(result["content"][0]["text"].as_str().unwrap().len() <= 16_384);
        assert!(
            result["structuredContent"].is_object(),
            "missing structured result"
        );
        result["structuredContent"].clone()
    }

    async fn tool_error(&mut self, name: &str, arguments: Value) {
        let response = self.call(name, arguments).await;
        assert!(
            response.get("error").is_none(),
            "expected application error, got {response}"
        );
        assert_eq!(
            response["result"]["isError"], true,
            "expected tool error: {response}"
        );
        let text = response["result"]["content"][0]["text"].as_str().unwrap();
        assert!(!text.is_empty() && text.len() <= 2048);
    }

    async fn finish(mut self) {
        // Stdio EOF, rather than a made-up shutdown method, ends MCP's lifecycle.
        drop(self.stdin.take());
        let output = timeout(DEADLINE, self.child.wait_with_output())
            .await
            .expect("server did not exit after stdin EOF")
            .unwrap();
        assert!(
            output.status.success(),
            "server exit failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stderr.is_empty(),
            "unexpected diagnostics: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let mut trailing = String::new();
        self.stdout.read_to_string(&mut trailing).await.unwrap();
        assert!(trailing.is_empty(), "unsolicited stdout: {trailing}");
    }
}

#[tokio::test]
async fn tools_have_typed_schemas_and_sdk_protocol_errors() {
    let fixture = Fixture::new();
    let mut server = Server::start(&fixture).await;
    let response = server.request("tools/list", json!({})).await;
    let tools = response["result"]["tools"].as_array().unwrap();
    let names: BTreeSet<_> = tools
        .iter()
        .map(|tool| tool["name"].as_str().unwrap())
        .collect();
    assert_eq!(
        names,
        BTreeSet::from([
            "dataset_open",
            "dataset_reopen",
            "dataset_preview",
            "dataset_query",
            "dataset_export",
            "dataset_read_artifact",
            "dataset_release",
            "storage_register",
            "storage_resolve",
            "prepared_fasta",
            "dataset_fasta_windows"
        ])
    );
    let prepared = tools
        .iter()
        .find(|tool| tool["name"] == "prepared_fasta")
        .unwrap();
    let required: BTreeSet<_> = prepared["inputSchema"]["required"]
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_str().unwrap())
        .collect();
    assert_eq!(
        required,
        BTreeSet::from(["reference", "snapshot_id", "source_path"])
    );
    let metrics = tools
        .iter()
        .find(|tool| tool["name"] == "dataset_fasta_windows")
        .unwrap();
    let required: BTreeSet<_> = metrics["inputSchema"]["required"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    assert_eq!(
        required,
        BTreeSet::from([
            "reference",
            "snapshot_id",
            "source_path",
            "recipe_id",
            "sequence_id",
            "window_size"
        ])
    );
    assert_eq!(metrics["annotations"]["readOnlyHint"], true);
    assert_eq!(prepared["annotations"]["readOnlyHint"], false);
    assert_eq!(prepared["annotations"]["destructiveHint"], false);
    for tool in tools {
        assert_eq!(tool["inputSchema"]["type"], "object");
        assert_eq!(tool["inputSchema"]["additionalProperties"], false);
        assert!(tool["inputSchema"]["properties"].is_object());
        assert_eq!(tool["annotations"]["openWorldHint"], false);
    }
    for (tool, arguments) in [
        ("dataset_open", json!({})),
        ("dataset_open", json!({"path": 12})),
        (
            "dataset_open",
            json!({"path": "genes.csv", "unexpected": true}),
        ),
        (
            "dataset_preview",
            json!({"dataset_id":"anything", "preview_rows":-1}),
        ),
        (
            "dataset_query",
            json!({"dataset_id":"anything", "filter":{"column":"score", "op":"execute", "value":"code"}}),
        ),
        ("storage_register", json!({})),
        (
            "storage_register",
            json!({"source_path":12, "requested_ref":"pdb:1ABC", "canonical_ref":"pdb:1ABC", "declaration":{"provider":"pdb","scope":"entry","representations":{"mmcif":["entry.cif"]}}}),
        ),
        (
            "storage_register",
            json!({"source_path":"pkg", "requested_ref":"pdb:1ABC", "canonical_ref":"pdb:1ABC", "declaration":{"provider":"invented"}}),
        ),
        ("storage_resolve", json!({"reference":"pdb:1ABC"})),
        (
            "storage_resolve",
            json!({"reference":"pdb:1ABC", "representation":"mmcif", "unexpected":true}),
        ),
        ("prepared_fasta", json!({})),
        (
            "prepared_fasta",
            json!({"reference":"refseq.gcf:GCF_000005845.2","snapshot_id":"sha256-anything"}),
        ),
        (
            "prepared_fasta",
            json!({"reference":"refseq.gcf:GCF_000005845.2","snapshot_id":42,"source_path":"genome.fna"}),
        ),
        (
            "prepared_fasta",
            json!({"reference":"refseq.gcf:GCF_000005845.2","snapshot_id":"sha256-anything","source_path":"genome.fna","unexpected":true}),
        ),
        ("dataset_fasta_windows", json!({})),
        (
            "dataset_fasta_windows",
            json!({"reference":"refseq.gcf:GCF_000005845.2","snapshot_id":"sha256-anything","source_path":"genome.fna","recipe_id":"sha256-anything","sequence_id":"chr", "window_size":-1}),
        ),
        (
            "dataset_fasta_windows",
            json!({"reference":"refseq.gcf:GCF_000005845.2","snapshot_id":"sha256-anything","source_path":"genome.fna","recipe_id":"sha256-anything","sequence_id":"chr", "window_size":1,"unexpected":true}),
        ),
        ("no_such_tool", json!({})),
    ] {
        let response = server.call(tool, arguments).await;
        assert_eq!(
            response["error"]["code"], -32602,
            "SDK invalid-params handling: {response}"
        );
    }
    // Bad requests must not kill an initialized session.
    let opened = server.successful("dataset_open", json!({"path":"genes.csv", "schema":{"score":"int64", "active":"boolean"}, "preview_rows":0})).await;
    assert_eq!(opened["row_count"], 5);
    server.finish().await;
}

#[tokio::test]
async fn full_data_query_export_retrieval_and_release_round_trip() {
    let fixture = Fixture::new();
    let mut server = Server::start(&fixture).await;
    let opened = server.successful("dataset_open", json!({"path":"genes.csv", "schema":{"score":"int64", "active":"boolean"}, "preview_rows":1})).await;
    let dataset_id = opened["dataset_id"].as_str().unwrap();
    assert!(dataset_id.starts_with("ds_"));
    assert_eq!(opened["row_count"], 5);
    assert_eq!(opened["engine"], "rust_polars");
    assert_eq!(
        opened["preview"]["rows"],
        json!([["BRCA1", 5, "low", true]])
    );
    assert_eq!(opened["preview"]["omitted_rows"], 4);
    assert_eq!(opened["preview"]["not_complete_data"], true);
    for key in ["identifier", "species", "reference", "coordinates", "units"] {
        assert!(
            opened["scientific_metadata"][key].is_null(),
            "must not infer {key}"
        );
    }
    assert_eq!(
        opened["source_sha256"],
        format!("{:x}", Sha256::digest(CSV.as_bytes()))
    );

    let preview = server
        .successful(
            "dataset_preview",
            json!({"dataset_id":dataset_id, "preview_rows":2}),
        )
        .await;
    assert_eq!(preview["dataset_id"], dataset_id);
    assert_eq!(preview["preview"]["returned_rows"], 2);
    let queried = server.successful("dataset_query", json!({
        "dataset_id":dataset_id, "filter":{"column":"score", "op":"gt", "value":10},
        "sort":{"column":"score", "descending":true}, "select":["gene","score"], "preview_rows":1
    })).await;
    let derived_id = queried["dataset_id"].as_str().unwrap();
    assert_ne!(derived_id, dataset_id);
    assert_eq!(
        queried["row_count"], 3,
        "query must use full data, not the one-row preview"
    );
    assert_eq!(queried["preview"]["rows"], json!([["TP53", 30]]));
    let nulls = server
        .successful(
            "dataset_query",
            json!({"dataset_id":dataset_id, "filter":{"column":"score", "op":"is_null"}}),
        )
        .await;
    assert_eq!(nulls["row_count"], 1);
    assert_eq!(nulls["preview"]["rows"][0][0], "KRAS");
    assert!(nulls["preview"]["rows"][0][1].is_null());

    let artifact = server
        .successful("dataset_export", json!({"dataset_id":derived_id}))
        .await;
    let artifact_id = artifact["artifact_id"].as_str().unwrap();
    assert_eq!(artifact["format"], "arrow_ipc");
    assert_eq!(artifact["record"]["row_count"], 3);
    assert_eq!(
        artifact["record"]["provenance"]["source_sha256"],
        opened["source_sha256"]
    );
    assert_eq!(
        artifact["record"]["provenance"]["operations"][0]["parent_dataset_id"],
        dataset_id
    );
    let mut bytes = Vec::new();
    let mut offset = 0_u64;
    loop {
        let chunk = server
            .successful(
                "dataset_read_artifact",
                json!({"artifact_id":artifact_id, "offset":offset, "max_bytes":97}),
            )
            .await;
        assert_eq!(chunk["offset"], offset);
        assert_eq!(chunk["encoding"], "base64");
        assert_eq!(chunk["sha256"], artifact["sha256"]);
        let decoded = STANDARD.decode(chunk["data"].as_str().unwrap()).unwrap();
        assert!(decoded.len() <= 97);
        bytes.extend(decoded);
        let next = chunk["next_offset"].as_u64().unwrap();
        assert_eq!(next as usize, bytes.len());
        if chunk["eof"] == true {
            break;
        }
        assert!(next > offset, "artifact retrieval must make progress");
        offset = next;
    }
    assert_eq!(bytes.len() as u64, artifact["bytes"].as_u64().unwrap());
    assert_eq!(format!("{:x}", Sha256::digest(&bytes)), artifact["sha256"]);
    let frame = IpcReader::new(Cursor::new(bytes.clone())).finish().unwrap();
    assert_eq!(frame.shape(), (3, 2));
    assert_eq!(
        frame
            .column("gene")
            .unwrap()
            .str()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>(),
        vec!["TP53", "ALK", "EGFR"]
    );
    assert_eq!(
        frame
            .column("score")
            .unwrap()
            .i64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>(),
        vec![30, 30, 20]
    );
    let artifact_path = PathBuf::from(artifact["execution_host_path"].as_str().unwrap());
    assert!(artifact_path.starts_with(&fixture.output));
    assert_eq!(fs::read(&artifact_path).unwrap(), bytes);
    let record: Value =
        serde_json::from_slice(&fs::read(artifact["record_path"].as_str().unwrap()).unwrap())
            .unwrap();
    assert_eq!(record, artifact["record"]);

    let released = server
        .successful("dataset_release", json!({"dataset_id":derived_id}))
        .await;
    assert_eq!(released["released"], derived_id);
    server
        .tool_error("dataset_preview", json!({"dataset_id":derived_id}))
        .await;
    server
        .tool_error("dataset_release", json!({"dataset_id":derived_id}))
        .await;
    assert!(
        artifact_path.is_file(),
        "release must preserve exported files"
    );
    server
        .successful(
            "dataset_read_artifact",
            json!({"artifact_id":artifact_id, "max_bytes":1}),
        )
        .await;
    server.finish().await;

    // Handles are capabilities within a single process, not persistent IDs.
    let mut restarted = Server::start(&fixture).await;
    restarted
        .tool_error("dataset_preview", json!({"dataset_id":dataset_id}))
        .await;
    restarted
        .tool_error(
            "dataset_read_artifact",
            json!({"artifact_id":artifact_id, "max_bytes":1}),
        )
        .await;
    restarted.finish().await;
}

#[tokio::test]
async fn validation_errors_are_bounded_tool_results_and_keep_session_alive() {
    let fixture = Fixture::new();
    let mut server = Server::start(&fixture).await;
    for path in [
        "../genes.csv",
        "/etc/passwd",
        "https://example.org/genes.csv",
        "missing.csv",
    ] {
        server
            .tool_error("dataset_open", json!({"path":path}))
            .await;
    }
    server
        .tool_error(
            "dataset_open",
            json!({"path":"genes.csv", "preview_rows":51}),
        )
        .await;
    server
        .tool_error(
            "dataset_open",
            json!({"path":"genes.csv", "metadata":{"identifier":"unsupported:example"}}),
        )
        .await;
    server
        .tool_error("dataset_preview", json!({"dataset_id":"unknown"}))
        .await;
    server
        .tool_error(
            "dataset_read_artifact",
            json!({"artifact_id":"unknown", "max_bytes":1}),
        )
        .await;
    let opened = server.successful("dataset_open", json!({"path":"genes.csv", "schema":{"score":"int64", "active":"boolean"}, "preview_rows":0})).await;
    let id = &opened["dataset_id"];
    for arguments in [
        json!({"dataset_id":id, "filter":{"column":"score", "op":"gt", "value":"10"}}),
        json!({"dataset_id":id, "filter":{"column":"active", "op":"gt", "value":true}}),
        json!({"dataset_id":id, "filter":{"column":"score", "op":"is_null", "value":0}}),
        json!({"dataset_id":id, "select":["gene","gene"]}),
        json!({"dataset_id":id, "select":["missing"]}),
        json!({"dataset_id":id, "preview_rows":51}),
    ] {
        server.tool_error("dataset_query", arguments).await;
    }
    let artifact = server
        .successful("dataset_export", json!({"dataset_id":id}))
        .await;
    for arguments in [
        json!({"artifact_id":artifact["artifact_id"], "max_bytes":0}),
        json!({"artifact_id":artifact["artifact_id"], "max_bytes":49153}),
        json!({"artifact_id":artifact["artifact_id"], "max_bytes":1, "offset":artifact["bytes"].as_u64().unwrap()+1}),
    ] {
        server.tool_error("dataset_read_artifact", arguments).await;
    }
    fs::write(
        artifact["execution_host_path"].as_str().unwrap(),
        b"changed",
    )
    .unwrap();
    server
        .tool_error(
            "dataset_read_artifact",
            json!({"artifact_id":artifact["artifact_id"], "max_bytes":1}),
        )
        .await;
    server
        .successful("dataset_preview", json!({"dataset_id":id}))
        .await;
    server.finish().await;
}

#[tokio::test]
async fn default_string_columns_preserve_identifiers_and_large_integer_lexemes() {
    let fixture = Fixture::new();
    fs::write(
        fixture.data.join("identifiers.csv"),
        "accession,wide\n00123,9007199254740993\n00001,9007199254740994\n",
    )
    .unwrap();
    let mut server = Server::start(&fixture).await;
    let opened = server
        .successful("dataset_open", json!({"path":"identifiers.csv"}))
        .await;
    assert_eq!(
        opened["preview"]["rows"],
        json!([["00123", "9007199254740993"], ["00001", "9007199254740994"]])
    );
    let typed = server
        .successful(
            "dataset_open",
            json!({"path":"identifiers.csv", "schema":{"wide":"int64"}}),
        )
        .await;
    assert_eq!(
        typed["preview"]["rows"][0],
        json!(["00123", 9007199254740993_i64])
    );
    let selected = server.successful("dataset_query", json!({"dataset_id":typed["dataset_id"], "filter":{"column":"wide", "op":"eq", "value":9007199254740993_i64}})).await;
    assert_eq!(selected["row_count"], 1);
    assert_eq!(selected["preview"]["rows"][0][0], "00123");
    server
        .tool_error(
            "dataset_open",
            json!({"path":"identifiers.csv", "schema":{"missing":"int64"}}),
        )
        .await;
    server.finish().await;
}

#[tokio::test]
async fn preview_and_maximum_artifact_chunks_respect_wire_budgets() {
    let fixture = Fixture::new();
    let header = (0..32)
        .map(|index| format!("column{index}"))
        .collect::<Vec<_>>()
        .join(",");
    let row = vec!["x".repeat(300); 32].join(",");
    fs::write(
        fixture.data.join("large.csv"),
        format!("{header}\n{}\n", vec![row; 8].join("\n")),
    )
    .unwrap();
    let mut server = Server::start(&fixture).await;
    let opened = server
        .successful(
            "dataset_open",
            json!({"path":"large.csv", "preview_rows":50}),
        )
        .await;
    assert_eq!(opened["row_count"], 8);
    assert!(opened["preview"]["returned_rows"].as_u64().unwrap() < 8);
    assert!(opened["preview"]["truncated_cells"].as_u64().unwrap() > 0);
    assert_eq!(opened["preview"]["not_complete_data"], true);
    for row in opened["preview"]["rows"].as_array().unwrap() {
        for cell in row.as_array().unwrap() {
            assert!(cell.as_str().unwrap().len() <= 256);
        }
    }
    let artifact = server
        .successful("dataset_export", json!({"dataset_id":opened["dataset_id"]}))
        .await;
    let chunk = server
        .successful(
            "dataset_read_artifact",
            json!({"artifact_id":artifact["artifact_id"], "max_bytes":49152}),
        )
        .await;
    assert_eq!(
        STANDARD
            .decode(chunk["data"].as_str().unwrap())
            .unwrap()
            .len(),
        49152
    );
    assert_eq!(chunk["eof"], false);
    server.finish().await;
}

#[cfg(unix)]
#[tokio::test]
async fn symlink_escape_is_rejected_over_real_stdio() {
    use std::os::unix::fs::symlink;
    let fixture = Fixture::new();
    fs::write(fixture.root.join("outside.csv"), CSV).unwrap();
    symlink(
        fixture.root.join("outside.csv"),
        fixture.data.join("escape.csv"),
    )
    .unwrap();
    let mut server = Server::start(&fixture).await;
    server
        .tool_error("dataset_open", json!({"path":"escape.csv"}))
        .await;
    server.finish().await;
}

#[tokio::test]
async fn command_help_and_errors_never_write_to_stdout() {
    for (args, success) in [
        (vec!["--help"], true),
        (vec!["mcp-native", "--help"], true),
        (vec!["--version"], true),
        (vec![], false),
        (vec!["unknown"], false),
        (vec!["mcp-native", "--unknown"], false),
        (vec!["mcp-native", "--data-root", "missing"], false),
    ] {
        let output = timeout(DEADLINE, Command::new(binary()).args(args).output())
            .await
            .unwrap()
            .unwrap();
        assert_eq!(output.status.success(), success);
        assert!(
            output.stdout.is_empty(),
            "CLI diagnostics must never corrupt MCP stdout"
        );
        assert!(!output.stderr.is_empty());
    }
    let fixture = Fixture::new();
    let output = Command::new(binary())
        .arg("mcp-native")
        .arg("--data-root")
        .arg(fixture.root.join("missing"))
        .arg("--output-root")
        .arg(&fixture.output)
        .output()
        .await
        .unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(!output.stderr.is_empty());
}

async fn downloaded_frame(server: &mut Server, exported: &Value) -> DataFrame {
    let mut bytes = Vec::new();
    loop {
        let chunk = server
            .successful(
                "dataset_read_artifact",
                json!({
                    "artifact_id":exported["artifact_id"], "offset":bytes.len(), "max_bytes":997
                }),
            )
            .await;
        assert_eq!(chunk["offset"], bytes.len());
        assert_eq!(chunk["sha256"], exported["sha256"]);
        bytes.extend(STANDARD.decode(chunk["data"].as_str().unwrap()).unwrap());
        if chunk["eof"] == true {
            break;
        }
    }
    assert_eq!(format!("{:x}", Sha256::digest(&bytes)), exported["sha256"]);
    assert_eq!(bytes.len(), exported["bytes"].as_u64().unwrap() as usize);
    IpcReader::new(Cursor::new(bytes)).finish().unwrap()
}

#[tokio::test]
async fn separate_processes_reopen_preserve_complete_typed_data_and_metadata() {
    let fixture = Fixture::new();
    fs::write(fixture.data.join("typed.csv"), "sample,count,ratio,active\n00123,9007199254740993,1.25,true\n00002,,,false\n00002,,,false\n00003,20,2.5,\n00004,-9,-0.0,true\nskip,1,0.0,false\n").unwrap();
    let metadata = json!({"identifier":"refseq.gcf:GCF_000001405.40", "species":null, "reference":"caller reference", "coordinates":null, "units":"caller count"});
    let mut first = Server::start(&fixture).await;
    let opened = first.successful("dataset_open",json!({"path":"typed.csv","schema":{"count":"int64","ratio":"float64","active":"boolean"},"metadata":metadata,"preview_rows":1})).await;
    let derived = first.successful("dataset_query",json!({"dataset_id":opened["dataset_id"],"filter":{"column":"sample","op":"lt","value":"01000"},"sort":{"column":"sample","descending":false},"select":["sample","count","ratio","active"],"preview_rows":1})).await;
    assert_eq!(derived["row_count"], 5);
    let exported = first
        .successful(
            "dataset_export",
            json!({"dataset_id":derived["dataset_id"]}),
        )
        .await;
    let original = downloaded_frame(&mut first, &exported).await;
    assert_eq!(original.height(), 5);
    assert_eq!(
        original
            .column("sample")
            .unwrap()
            .str()
            .unwrap()
            .into_iter()
            .collect::<Vec<_>>(),
        vec![
            Some("00002"),
            Some("00002"),
            Some("00003"),
            Some("00004"),
            Some("00123")
        ]
    );
    assert_eq!(
        original.column("count").unwrap().i64().unwrap().get(4),
        Some(9007199254740993_i64)
    );
    assert_eq!(original.column("count").unwrap().null_count(), 2);
    first.finish().await;
    let second_output = fixture.root.join("second-output");
    fs::create_dir(&second_output).unwrap();
    let record_name = Path::new(exported["record_path"].as_str().unwrap())
        .file_name()
        .unwrap()
        .to_str()
        .unwrap();
    let mut second = Server::start_roots(&fixture.output, &second_output).await;
    second
        .tool_error(
            "dataset_preview",
            json!({"dataset_id":derived["dataset_id"]}),
        )
        .await;
    second
        .tool_error(
            "dataset_read_artifact",
            json!({"artifact_id":exported["artifact_id"],"offset":0,"max_bytes":10}),
        )
        .await;
    let reopened = second
        .successful(
            "dataset_reopen",
            json!({"record_path":record_name,"preview_rows":1}),
        )
        .await;
    assert_ne!(reopened["dataset_id"], derived["dataset_id"]);
    assert_eq!(reopened["schema"], derived["schema"]);
    assert_eq!(reopened["row_count"], 5);
    assert_eq!(reopened["preview"]["omitted_rows"], 4);
    assert_eq!(reopened["scientific_metadata"], metadata);
    assert_eq!(reopened["source_sha256"], opened["source_sha256"]);
    assert_eq!(
        reopened["reopen_verification"]["artifact_sha256"],
        exported["sha256"]
    );
    assert_eq!(
        reopened["reopen_verification"]["record_sha256"],
        format!(
            "{:x}",
            Sha256::digest(fs::read(fixture.output.join(record_name)).unwrap())
        )
    );
    assert_eq!(
        reopened["reopen_verification"]["authenticity"],
        "not_established"
    );
    assert_eq!(
        reopened["reopen_verification"]["original_provenance"],
        "recorded_claims_not_independently_verified"
    );
    let reexported = second
        .successful(
            "dataset_export",
            json!({"dataset_id":reopened["dataset_id"]}),
        )
        .await;
    let all_rows = downloaded_frame(&mut second, &reexported).await;
    assert_eq!(original.schema(), all_rows.schema());
    assert!(original
        .column("ratio")
        .unwrap()
        .f64()
        .unwrap()
        .get(3)
        .unwrap()
        .is_sign_negative());
    assert_eq!(
        original
            .column("ratio")
            .unwrap()
            .f64()
            .unwrap()
            .into_iter()
            .map(|value| value.map(f64::to_bits))
            .collect::<Vec<_>>(),
        all_rows
            .column("ratio")
            .unwrap()
            .f64()
            .unwrap()
            .into_iter()
            .map(|value| value.map(f64::to_bits))
            .collect::<Vec<_>>()
    );
    assert!(
        original.equals_missing(&all_rows),
        "every value, null, dtype and row must survive reopening"
    );
    assert_eq!(reexported["record"]["scientific_metadata"], metadata);
    assert_eq!(
        reexported["record"]["reopen_verification"],
        reopened["reopen_verification"]
    );
    assert_eq!(
        reexported["record"]["provenance"],
        exported["record"]["provenance"]
    );
    let filtered = second.successful("dataset_query",json!({"dataset_id":reopened["dataset_id"],"filter":{"column":"count","op":"eq","value":9007199254740993_i64},"preview_rows":1})).await;
    assert_eq!(
        filtered["preview"]["rows"],
        json!([["00123", 9007199254740993_i64, 1.25, true]])
    );
    let filtered_export = second
        .successful(
            "dataset_export",
            json!({"dataset_id":filtered["dataset_id"]}),
        )
        .await;
    assert_eq!(
        downloaded_frame(&mut second, &filtered_export)
            .await
            .height(),
        1
    );
    second.finish().await;
    let third_output = fixture.root.join("third-output");
    fs::create_dir(&third_output).unwrap();
    let mut third = Server::start_roots(&second_output, &third_output).await;
    let record_name = Path::new(reexported["record_path"].as_str().unwrap())
        .file_name()
        .unwrap()
        .to_str()
        .unwrap();
    let reopened_again = third
        .successful("dataset_reopen", json!({"record_path":record_name}))
        .await;
    assert_eq!(reopened_again["row_count"], 5);
    assert_eq!(reopened_again["scientific_metadata"], metadata);
    third.finish().await;
}

#[tokio::test]
async fn reopen_requires_matching_record_and_rejects_path_and_schema_errors() {
    let fixture = Fixture::new();
    let mut first = Server::start(&fixture).await;
    let opened = first
        .successful(
            "dataset_open",
            json!({"path":"genes.csv","schema":{"score":"int64","active":"boolean"}}),
        )
        .await;
    let exported = first
        .successful("dataset_export", json!({"dataset_id":opened["dataset_id"]}))
        .await;
    first.finish().await;
    let record_path = PathBuf::from(exported["record_path"].as_str().unwrap());
    let record_name = record_path.file_name().unwrap().to_str().unwrap();
    let record_bytes = fs::read(&record_path).unwrap();
    let base_record: Value = serde_json::from_slice(&record_bytes).unwrap();
    let ipc_path = PathBuf::from(exported["execution_host_path"].as_str().unwrap());
    let ipc_bytes = fs::read(&ipc_path).unwrap();
    let mut second = Server::start_roots(&fixture.output, &fixture.data).await;
    for path in [
        "missing.json",
        "../missing.json",
        record_path.to_str().unwrap(),
    ] {
        second
            .tool_error("dataset_reopen", json!({"record_path":path}))
            .await;
    }
    second
        .tool_error(
            "dataset_reopen",
            json!({"record_path":record_name,"preview_rows":51}),
        )
        .await;
    let mut corruptions = Vec::new();
    let mut value = base_record.clone();
    value["sha256"] = json!("0".repeat(64));
    corruptions.push(value);
    let mut value = base_record.clone();
    value["row_count"] = json!(999);
    corruptions.push(value);
    let mut value = base_record.clone();
    value["schema"][0]["dtype"] = json!("i64");
    corruptions.push(value);
    let mut value = base_record.clone();
    value["schema"][0]["name"] = json!("different");
    corruptions.push(value);
    let mut value = base_record.clone();
    value["file"] = json!("../genes.csv");
    corruptions.push(value);
    let mut value = base_record.clone();
    value["record_version"] = json!(999);
    corruptions.push(value);
    let mut value = base_record.clone();
    value["scientific_metadata"]["identifier"] = json!("ds_not_biological");
    corruptions.push(value);
    for record in corruptions {
        fs::write(&record_path, serde_json::to_vec(&record).unwrap()).unwrap();
        second
            .tool_error("dataset_reopen", json!({"record_path":record_name}))
            .await;
    }
    fs::write(&record_path, &record_bytes).unwrap();
    fs::remove_file(&ipc_path).unwrap();
    second
        .tool_error("dataset_reopen", json!({"record_path":record_name}))
        .await;
    fs::write(&ipc_path, b"tampered artifact").unwrap();
    second
        .tool_error("dataset_reopen", json!({"record_path":record_name}))
        .await;
    fs::write(&ipc_path, &ipc_bytes).unwrap();
    fs::remove_file(&record_path).unwrap();
    second
        .tool_error("dataset_reopen", json!({"record_path":record_name}))
        .await;
    fs::write(&record_path, &record_bytes).unwrap();
    let reopened = second
        .successful("dataset_reopen", json!({"record_path":record_name}))
        .await;
    assert_eq!(reopened["row_count"], 5);
    assert!(reopened["scientific_metadata"]["species"].is_null());
    second.finish().await;
}

// The identity and representation are intentionally caller declarations. This
// tiny synthetic file is not claimed to be an authenticated PDB provider record.
const DECLARED_MMCIF: &str = "data_synthetic\n_entry.id 1ABC\n#\n";

fn declared_pdb_request(source_path: &str) -> Value {
    json!({
        "source_path":source_path,
        "requested_ref":"pdb:1abc",
        "canonical_ref":"pdb:1ABC",
        "declaration":{
            "provider":"pdb", "scope":"entry",
            "representations":{"mmcif":["entry.cif"]}
        }
    })
}

fn prepare_native_fixture(fixture: &Fixture) -> PathBuf {
    fs::create_dir(fixture.data.join("native")).unwrap();
    fs::write(fixture.data.join("native/entry.cif"), DECLARED_MMCIF).unwrap();
    let store = fixture.root.join("native-store");
    fs::create_dir(&store).unwrap();
    store
}

#[tokio::test]
async fn storage_tools_report_not_configured_without_affecting_datasets() {
    let fixture = Fixture::new();
    let mut server = Server::start(&fixture).await;
    for (tool, request) in [
        ("storage_register", declared_pdb_request("native")),
        (
            "storage_resolve",
            json!({"reference":"pdb:1ABC","representation":"mmcif"}),
        ),
        (
            "prepared_fasta",
            prepared_request(&format!("sha256-{}", "a".repeat(64))),
        ),
    ] {
        let response = server.call(tool, request).await;
        assert_eq!(response["result"]["isError"], true);
        assert!(response["result"]["content"][0]["text"]
            .as_str()
            .unwrap()
            .contains("not configured"));
        assert!(response["result"]["content"][0]["text"]
            .as_str()
            .unwrap()
            .contains("--store-root"));
    }
    let opened = server
        .successful("dataset_open", json!({"path":"genes.csv"}))
        .await;
    assert_eq!(opened["row_count"], 5);
    server.finish().await;
}

#[tokio::test]
async fn storage_register_resolve_survives_process_restart_and_store_move() {
    let fixture = Fixture::new();
    let store = prepare_native_fixture(&fixture);
    let mut first = Server::start_with_storage(&fixture.data, &fixture.output, Some(&store)).await;
    let registered = first
        .successful("storage_register", declared_pdb_request("native"))
        .await;
    assert_eq!(registered["reused"], false);
    assert_eq!(registered["requested_ref"], "pdb:1abc");
    assert_eq!(registered["snapshot"]["canonical_ref"], "pdb:1ABC");
    assert_eq!(registered["snapshot"]["scope"], "entry");
    assert_eq!(registered["snapshot"]["file_count"], 1);
    assert!(registered.get("inventory").is_none());
    assert!(registered["snapshot"].get("inventory").is_none());
    let snapshot_id = registered["snapshot"]["snapshot_id"].clone();
    let request =
        json!({"reference":"pdb:1ABC","representation":"mmcif","snapshot_id":snapshot_id});
    let resolved = first.successful("storage_resolve", request.clone()).await;
    assert_eq!(resolved["status"], "ready");
    assert_eq!(resolved["paths"].as_array().unwrap().len(), 1);
    let relative = resolved["paths"][0]["relative_path"].as_str().unwrap();
    assert!(!Path::new(relative).is_absolute());
    assert_eq!(
        fs::read_to_string(store.join(relative)).unwrap(),
        DECLARED_MMCIF
    );
    assert_eq!(
        fs::read_to_string(
            resolved["paths"][0]["execution_host_path"]
                .as_str()
                .unwrap()
        )
        .unwrap(),
        DECLARED_MMCIF
    );
    assert!(resolved.get("inventory").is_none());
    assert!(resolved.get("data").is_none());
    let reused = first
        .successful("storage_register", declared_pdb_request("native"))
        .await;
    assert_eq!(reused["reused"], true);
    assert_eq!(reused["snapshot"]["snapshot_id"], snapshot_id);
    let opened = first
        .successful("dataset_open", json!({"path":"genes.csv"}))
        .await;
    assert_eq!(opened["row_count"], 5);
    first.finish().await;

    fs::remove_dir_all(&fixture.data).unwrap();
    let moved = fixture.root.join("moved-store");
    fs::rename(&store, &moved).unwrap();
    let empty_data = fixture.root.join("unrelated-empty-data");
    fs::create_dir(&empty_data).unwrap();
    let mut second = Server::start_with_storage(&empty_data, &fixture.output, Some(&moved)).await;
    let reopened = second.successful("storage_resolve", request).await;
    assert_eq!(reopened["status"], "ready");
    assert_eq!(reopened["snapshot"]["snapshot_id"], snapshot_id);
    assert_eq!(reopened["paths"][0]["relative_path"], relative);
    let host_path = Path::new(
        reopened["paths"][0]["execution_host_path"]
            .as_str()
            .unwrap(),
    );
    assert!(host_path.starts_with(&moved));
    assert_eq!(fs::read_to_string(host_path).unwrap(), DECLARED_MMCIF);
    let unknown = second
        .successful(
            "storage_resolve",
            json!({"reference":"pdb:1ABC", "representation":"not_supplied"}),
        )
        .await;
    assert_eq!(unknown["status"], "unavailable");
    assert_eq!(unknown["available_representations"], json!(["mmcif"]));
    let absent = second
        .successful(
            "storage_resolve",
            json!({"reference":"pdb:9XYZ", "representation":"mmcif"}),
        )
        .await;
    assert_eq!(absent["status"], "miss");
    // Full integrity checking is observable as a structured status, not raw bytes.
    fs::write(host_path, b"changed native bytes\n").unwrap();
    let corrupt = second
        .successful(
            "storage_resolve",
            json!({"reference":"pdb:1ABC","representation":"mmcif","snapshot_id":snapshot_id}),
        )
        .await;
    assert_eq!(corrupt["status"], "corrupt");
    second.finish().await;
}

#[tokio::test]
async fn storage_invalid_paths_ids_and_scopes_are_bounded_errors() {
    let fixture = Fixture::new();
    let store = prepare_native_fixture(&fixture);
    let mut server = Server::start_with_storage(&fixture.data, &fixture.output, Some(&store)).await;
    for path in [
        "../native",
        "/etc",
        "missing",
        "https://example.org/package",
    ] {
        server
            .tool_error("storage_register", declared_pdb_request(path))
            .await;
    }
    for changes in [
        json!({"requested_ref":"pdb:2ABC"}),
        json!({"canonical_ref":"ds_session_handle"}),
        json!({"declaration":{"provider":"pdb","scope":"assembly:0","representations":{"mmcif":["entry.cif"]}}}),
        json!({"declaration":{"provider":"pdb","scope":"entry","representations":{"mmcif":["../entry.cif"]}}}),
        json!({"declaration":{"provider":"pdb","scope":"entry","representations":{"mmcif":["missing.cif"]}}}),
    ] {
        let mut request = declared_pdb_request("native");
        for (key, value) in changes.as_object().unwrap() {
            request[key] = value.clone();
        }
        server.tool_error("storage_register", request).await;
    }
    for request in [
        json!({"reference":"ds_session_handle", "representation":"mmcif"}),
        json!({"reference":"pdb:1ABC", "representation":"../mmcif"}),
        json!({"reference":"pdb:1ABC", "representation":"mmcif", "snapshot_id":"../snapshot"}),
        json!({"reference":"pdb:1ABC", "representation":"mmcif", "scope":"assembly:0"}),
    ] {
        server.tool_error("storage_resolve", request).await;
    }
    let registered = server
        .successful("storage_register", declared_pdb_request("native"))
        .await;
    let other_scope = server.successful("storage_resolve", json!({"reference":"pdb:1ABC", "representation":"mmcif", "scope":"assembly:1", "snapshot_id":registered["snapshot"]["snapshot_id"]})).await;
    assert_eq!(
        other_scope["status"], "miss",
        "entry does not imply assembly scope"
    );
    server.finish().await;
}

#[tokio::test]
async fn storage_ambiguity_requires_explicit_snapshot_selection() {
    let fixture = Fixture::new();
    let store = prepare_native_fixture(&fixture);
    let mut server = Server::start_with_storage(&fixture.data, &fixture.output, Some(&store)).await;
    let first = server
        .successful("storage_register", declared_pdb_request("native"))
        .await;
    fs::write(
        fixture.data.join("native/entry.cif"),
        "data_other_synthetic\n_entry.id 1ABC\n#\n",
    )
    .unwrap();
    let second = server
        .successful("storage_register", declared_pdb_request("native"))
        .await;
    assert_ne!(
        first["snapshot"]["snapshot_id"],
        second["snapshot"]["snapshot_id"]
    );
    let ambiguous = server
        .successful(
            "storage_resolve",
            json!({"reference":"pdb:1ABC", "representation":"mmcif"}),
        )
        .await;
    assert_eq!(ambiguous["status"], "ambiguous");
    assert_eq!(ambiguous["candidates"].as_array().unwrap().len(), 2);
    let selected = server.successful("storage_resolve", json!({"reference":"pdb:1ABC", "representation":"mmcif", "snapshot_id":first["snapshot"]["snapshot_id"]})).await;
    assert_eq!(selected["status"], "ready");
    assert_eq!(
        fs::read_to_string(
            selected["paths"][0]["execution_host_path"]
                .as_str()
                .unwrap()
        )
        .unwrap(),
        DECLARED_MMCIF
    );
    server.finish().await;
}

async fn cli_storage(
    fixture: &Fixture,
    store: &Path,
    action: &str,
    request_file: &Path,
) -> std::process::Output {
    let mut command = Command::new(binary());
    command
        .args(["storage", action])
        .arg("--store-root")
        .arg(store)
        .arg("--request-file")
        .arg(request_file);
    if action == "register" {
        command.arg("--source-root").arg(&fixture.data);
    }
    timeout(
        DEADLINE,
        command.env("PATH", "").env_remove("PYTHONPATH").output(),
    )
    .await
    .unwrap()
    .unwrap()
}

fn successful_cli_json(output: &std::process::Output) -> Value {
    assert!(
        output.status.success(),
        "CLI failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    assert!(output.stdout.len() <= 64 * 1024 + 1);
    let result: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(result.is_object());
    result
}

#[tokio::test]
async fn one_shot_storage_cli_has_bounded_strict_json_and_no_source_for_resolution() {
    let fixture = Fixture::new();
    let store = prepare_native_fixture(&fixture);
    let request_file = fixture.root.join("request.json");
    fs::write(
        &request_file,
        serde_json::to_vec(&declared_pdb_request("native")).unwrap(),
    )
    .unwrap();
    let registered =
        successful_cli_json(&cli_storage(&fixture, &store, "register", &request_file).await);
    assert_eq!(registered["reused"], false);
    fs::remove_dir_all(&fixture.data).unwrap();
    let moved = fixture.root.join("moved-cli-store");
    fs::rename(&store, &moved).unwrap();
    let request = json!({"reference":"pdb:1ABC","representation":"mmcif"});
    fs::write(&request_file, serde_json::to_vec(&request).unwrap()).unwrap();
    let resolved =
        successful_cli_json(&cli_storage(&fixture, &moved, "resolve", &request_file).await);
    assert_eq!(resolved["status"], "ready");
    assert_eq!(
        resolved["snapshot"]["snapshot_id"],
        registered["snapshot"]["snapshot_id"]
    );
    assert_eq!(
        fs::read_to_string(
            resolved["paths"][0]["execution_host_path"]
                .as_str()
                .unwrap()
        )
        .unwrap(),
        DECLARED_MMCIF
    );
    for (reference, representation, status) in [
        ("pdb:9XYZ", "mmcif", "miss"),
        ("pdb:1ABC", "unknown_representation", "unavailable"),
    ] {
        fs::write(
            &request_file,
            serde_json::to_vec(&json!({"reference":reference,"representation":representation}))
                .unwrap(),
        )
        .unwrap();
        let result =
            successful_cli_json(&cli_storage(&fixture, &moved, "resolve", &request_file).await);
        assert_eq!(result["status"], status);
    }
    for invalid in [
        "{}".to_owned(),
        "{broken JSON".to_owned(),
        r#"{"reference":"pdb:1ABC","reference":"pdb:2ABC","representation":"mmcif"}"#.to_owned(),
        r#"{"reference":"pdb:1ABC","representation":"mmcif","extra":true}"#.to_owned(),
        " ".repeat(64 * 1024 + 1),
    ] {
        fs::write(&request_file, invalid).unwrap();
        let output = cli_storage(&fixture, &moved, "resolve", &request_file).await;
        assert!(!output.status.success());
        assert!(
            output.stdout.is_empty(),
            "errors must not look like JSON results"
        );
        assert!(!output.stderr.is_empty() && output.stderr.len() <= 2060);
    }
    fs::remove_file(&request_file).unwrap();
    let output = cli_storage(&fixture, &moved, "resolve", &request_file).await;
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
}

#[cfg(unix)]
#[tokio::test]
async fn storage_cli_rejects_fifo_before_opening_it() {
    let fixture = Fixture::new();
    let store = prepare_native_fixture(&fixture);
    let request_file = fixture.root.join("request.fifo");
    // Test-only POSIX fixture creation. The actual biov subprocess still runs
    // with empty PATH and cannot depend on this or another helper executable.
    assert!(std::process::Command::new("mkfifo")
        .arg(&request_file)
        .status()
        .unwrap()
        .success());
    let output = cli_storage(&fixture, &store, "resolve", &request_file).await;
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("request file must be a regular file"));
}

#[tokio::test]
async fn csv_record_semantics_match_complete_retrieved_exports() {
    let fixture = Fixture::new();
    let cases = [
        ("x\ra\rb\r", json!([["a"], ["b"]])),
        ("x\r\na\r\nb\r\n", json!([["a"], ["b"]])),
        ("x,y\na,b\n\nc,d\n", json!([["a", "b"], ["c", "d"]])),
        ("x\n\"\"\n\n", json!([[""]])),
        (
            "\r\nx,y\r\"a\n\nb\",\"c\rd\"\r\n\"a\r\n\r\nb\",\"say \"\"hi\"\"\"\n,\"\"\r",
            json!([["a\n\nb", "c\rd"], ["a\r\n\r\nb", "say \"hi\""], [null, ""]]),
        ),
    ];
    let mut server = Server::start(&fixture).await;
    for (index, (csv, expected)) in cases.into_iter().enumerate() {
        let path = format!("syntax-{index}.csv");
        fs::write(fixture.data.join(&path), csv).unwrap();
        // No preview: every value below must come from the full export.
        let opened = server
            .successful("dataset_open", json!({"path":path,"preview_rows":0}))
            .await;
        assert_eq!(opened["row_count"], expected.as_array().unwrap().len());
        assert_eq!(
            opened["source_sha256"],
            format!("{:x}", Sha256::digest(csv.as_bytes()))
        );
        let exported = server
            .successful("dataset_export", json!({"dataset_id":opened["dataset_id"]}))
            .await;
        assert_eq!(exported["record"]["row_count"], opened["row_count"]);
        let mut bytes = Vec::new();
        loop {
            let chunk = server.successful("dataset_read_artifact", json!({"artifact_id":exported["artifact_id"],"offset":bytes.len(),"max_bytes":113})).await;
            bytes.extend(STANDARD.decode(chunk["data"].as_str().unwrap()).unwrap());
            if chunk["eof"] == true {
                break;
            }
        }
        assert_eq!(format!("{:x}", Sha256::digest(&bytes)), exported["sha256"]);
        let frame = IpcReader::new(Cursor::new(bytes)).finish().unwrap();
        assert_eq!(frame.height(), expected.as_array().unwrap().len());
        for (row_index, row) in expected.as_array().unwrap().iter().enumerate() {
            assert_eq!(frame.width(), row.as_array().unwrap().len());
            for (column, value) in row.as_array().unwrap().iter().enumerate() {
                assert_eq!(
                    frame.get_columns()[column].str().unwrap().get(row_index),
                    value.as_str()
                );
            }
        }
    }
    for (index, csv) in ["x\n\"a\n", "x\n\"a\"junk\n", "x,y\ra\r", "x,y\r\na,b,c\r\n"]
        .into_iter()
        .enumerate()
    {
        let path = format!("invalid-{index}.csv");
        fs::write(fixture.data.join(&path), csv).unwrap();
        server
            .tool_error("dataset_open", json!({"path":path}))
            .await;
    }
    server.finish().await;
}

#[tokio::test]
async fn csv_blank_lines_and_retained_typed_payload_cannot_bypass_preflight() {
    let fixture = Fixture::new();
    let names: Vec<_> = (0..64).map(|i| format!("c{i}")).collect();
    let header = names.join(",");
    let row = format!("{}\n", vec!["1"; 64].join(","));
    fs::write(
        fixture.data.join("blank.csv"),
        format!("{header}\n{}{row}", "\n".repeat(62_000)),
    )
    .unwrap();
    let mut server = Server::start(&fixture).await;
    let blank = server
        .successful("dataset_open", json!({"path":"blank.csv","preview_rows":0}))
        .await;
    assert_eq!(blank["row_count"], 1);
    let exported = server
        .successful("dataset_export", json!({"dataset_id":blank["dataset_id"]}))
        .await;
    let chunk = server
        .successful(
            "dataset_read_artifact",
            json!({"artifact_id":exported["artifact_id"],"offset":0,"max_bytes":49152}),
        )
        .await;
    assert_eq!(chunk["eof"], true);
    let bytes = STANDARD.decode(chunk["data"].as_str().unwrap()).unwrap();
    let frame = IpcReader::new(Cursor::new(bytes)).finish().unwrap();
    assert_eq!(frame.shape(), (1, 64));
    assert!(frame
        .get_columns()
        .iter()
        .all(|column| column.str().unwrap().get(0) == Some("1")));
    server
        .successful("dataset_release", json!({"dataset_id":blank["dataset_id"]}))
        .await;

    let schema: serde_json::Map<String, Value> = names
        .iter()
        .map(|name| (name.clone(), json!("int64")))
        .collect();
    // Each table fits, but two exceed 64 MiB when numeric payload is included.
    fs::write(
        fixture.data.join("typed.csv"),
        format!("{header}\n{}", row.repeat(21_000)),
    )
    .unwrap();
    let first = server
        .successful(
            "dataset_open",
            json!({"path":"typed.csv","schema":schema,"preview_rows":0}),
        )
        .await;
    assert_eq!(first["row_count"], 21_000);
    let response = server
        .call(
            "dataset_open",
            json!({"path":"typed.csv","schema":schema,"preview_rows":0}),
        )
        .await;
    assert_eq!(response["result"]["isError"], true);
    assert!(response["result"]["content"][0]["text"]
        .as_str()
        .unwrap()
        .contains("allocation exceeds remaining retained dataset budget"));
    server
        .successful("dataset_release", json!({"dataset_id":first["dataset_id"]}))
        .await;
    let after_release = server
        .successful(
            "dataset_open",
            json!({"path":"typed.csv","schema":schema,"preview_rows":0}),
        )
        .await;
    assert_eq!(after_release["row_count"], 21_000);
    server.finish().await;
}

// Synthetic provider-shaped metadata is sufficient to exercise the native
// adapter. Supplied MD5 values establish consistency, never NCBI authenticity.
const PREPARED_REFERENCE: &str = "refseq.gcf:GCF_000005845.2";
const PREPARED_SOURCE: &str = "ncbi_dataset/data/GCF_000005845.2/genome.fna";
const PREPARED_FASTA: &str = ">chr1 synthetic fixture\nACGT\nAC\n>chr2\nttNN\n";

fn refseq_register_request() -> Value {
    json!({
        "source_path":"refseq",
        "requested_ref":"refseq.gcf:GCF_000005845",
        "canonical_ref":PREPARED_REFERENCE,
        "declaration":{"provider":"refseq"}
    })
}

fn prepared_request(snapshot_id: &str) -> Value {
    json!({"reference":PREPARED_REFERENCE,"snapshot_id":snapshot_id,"source_path":PREPARED_SOURCE})
}

fn prepare_refseq_fixture(fixture: &Fixture) -> PathBuf {
    let package = fixture.data.join("refseq");
    let data = package.join("ncbi_dataset/data");
    let assembly = data.join("GCF_000005845.2");
    fs::create_dir_all(&assembly).unwrap();
    fs::write(
        package.join("README.md"),
        "Synthetic native RefSeq-style fixture.\n",
    )
    .unwrap();
    fs::write(package.join(PREPARED_SOURCE), PREPARED_FASTA).unwrap();
    let report = b"{\"accession\":\"GCF_000005845.2\"}\n";
    fs::write(data.join("assembly_data_report.jsonl"), report).unwrap();
    let catalog = json!({"apiVersion":"V2","assemblies":[
        {"files":[{"filePath":"assembly_data_report.jsonl","fileType":"DATA_REPORT","uncompressedLengthBytes":report.len().to_string()}]},
        {"accession":"GCF_000005845.2","files":[
            {"filePath":"GCF_000005845.2/genome.fna","fileType":"GENOMIC_NUCLEOTIDE_FASTA","uncompressedLengthBytes":PREPARED_FASTA.len().to_string()}
        ]}
    ]});
    fs::write(
        data.join("dataset_catalog.json"),
        serde_json::to_vec_pretty(&catalog).unwrap(),
    )
    .unwrap();
    let mut md5 = String::new();
    for relative in [
        "ncbi_dataset/data/dataset_catalog.json",
        "ncbi_dataset/data/assembly_data_report.jsonl",
        PREPARED_SOURCE,
    ] {
        md5.push_str(&format!(
            "{:x}  {relative}\n",
            Md5::digest(fs::read(package.join(relative)).unwrap())
        ));
    }
    fs::write(package.join("md5sum.txt"), md5).unwrap();
    let store = fixture.root.join("native-store");
    fs::create_dir(&store).unwrap();
    store
}

fn file_hashes(root: &Path) -> std::collections::BTreeMap<PathBuf, Vec<u8>> {
    fn visit(
        root: &Path,
        directory: &Path,
        result: &mut std::collections::BTreeMap<PathBuf, Vec<u8>>,
    ) {
        for entry in fs::read_dir(directory).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                visit(root, &path, result);
            } else {
                result.insert(
                    path.strip_prefix(root).unwrap().to_path_buf(),
                    Sha256::digest(fs::read(path).unwrap()).to_vec(),
                );
            }
        }
    }
    let mut result = std::collections::BTreeMap::new();
    visit(root, root, &mut result);
    result
}

fn verify_prepared_paths(store: &Path, result: &Value) {
    assert_eq!(result["reference"], PREPARED_REFERENCE);
    assert_eq!(result["sequence_count"], 2);
    assert_eq!(result["total_bases"], 10);
    assert!(result["recipe_id"].is_string());
    assert!(result.get("sequences").is_none());
    assert!(result.get("inventory").is_none());
    assert!(result.get("data").is_none());
    for name in ["fasta", "fai", "dictionary", "provenance", "readme"] {
        let relative = result[name]["relative_path"].as_str().unwrap();
        let execution = Path::new(result[name]["execution_host_path"].as_str().unwrap());
        assert!(!Path::new(relative).is_absolute());
        assert_eq!(execution, store.join(relative));
        assert!(
            execution.is_file(),
            "missing {name} at {}",
            execution.display()
        );
    }
    assert_eq!(
        fs::read_to_string(result["fasta"]["execution_host_path"].as_str().unwrap()).unwrap(),
        PREPARED_FASTA
    );
    let fai = fs::read_to_string(result["fai"]["execution_host_path"].as_str().unwrap()).unwrap();
    assert_eq!(fai.lines().count(), 2);
    assert!(fai.starts_with("chr1\t6\t"));
    assert!(fai.lines().nth(1).unwrap().starts_with("chr2\t4\t"));
}

#[tokio::test]
async fn prepared_fasta_native_registration_reuse_and_moved_store_round_trip() {
    let fixture = Fixture::new();
    let store = prepare_refseq_fixture(&fixture);
    let native_original = file_hashes(&fixture.data.join("refseq"));
    let mut first = Server::start_with_storage(&fixture.data, &fixture.output, Some(&store)).await;
    let registered = first
        .successful("storage_register", refseq_register_request())
        .await;
    let snapshot = store.join(registered["snapshot"]["snapshot_path"].as_str().unwrap());
    let snapshot_before = file_hashes(&snapshot);
    let request = prepared_request(registered["snapshot"]["snapshot_id"].as_str().unwrap());
    let prepared = first.successful("prepared_fasta", request.clone()).await;
    assert_eq!(prepared["reused"], false);
    assert_eq!(
        prepared["snapshot_id"],
        registered["snapshot"]["snapshot_id"]
    );
    verify_prepared_paths(&store, &prepared);
    let repeated = first.successful("prepared_fasta", request.clone()).await;
    assert_eq!(repeated["reused"], true);
    assert_eq!(repeated["recipe_id"], prepared["recipe_id"]);
    for name in ["fasta", "fai", "dictionary", "provenance", "readme"] {
        assert_eq!(repeated[name], prepared[name]);
    }
    assert_eq!(
        file_hashes(&snapshot),
        snapshot_before,
        "native snapshot must remain immutable"
    );
    assert_eq!(
        file_hashes(&fixture.data.join("refseq")),
        native_original,
        "original package must remain immutable"
    );
    let opened = first
        .successful("dataset_open", json!({"path":"genes.csv"}))
        .await;
    assert_eq!(opened["row_count"], 5);
    first.finish().await;

    fs::remove_dir_all(fixture.data.join("refseq")).unwrap();
    let moved = fixture.root.join("moved-prepared-store");
    fs::rename(&store, &moved).unwrap();
    let mut second = Server::start_with_storage(&fixture.data, &fixture.output, Some(&moved)).await;
    let reused = second.successful("prepared_fasta", request).await;
    assert_eq!(reused["reused"], true);
    assert_eq!(reused["recipe_id"], prepared["recipe_id"]);
    verify_prepared_paths(&moved, &reused);
    for name in ["fasta", "fai", "dictionary", "provenance", "readme"] {
        assert_eq!(
            reused[name]["relative_path"],
            prepared[name]["relative_path"]
        );
    }
    second.finish().await;
}

#[tokio::test]
async fn prepared_fasta_invalid_reference_snapshot_and_paths_are_bounded_tool_errors() {
    let fixture = Fixture::new();
    let store = prepare_refseq_fixture(&fixture);
    let mut server = Server::start_with_storage(&fixture.data, &fixture.output, Some(&store)).await;
    let registered = server
        .successful("storage_register", refseq_register_request())
        .await;
    let request = prepared_request(registered["snapshot"]["snapshot_id"].as_str().unwrap());
    for (field, invalid) in [
        ("reference", "refseq.gcf:GCF_000005845".to_owned()),
        ("reference", "pdb:1ABC".to_owned()),
        ("reference", "ds_session_handle".to_owned()),
        ("reference", "refseq.gcf:GCF_000005845.3".to_owned()),
        ("snapshot_id", "../snapshot".to_owned()),
        ("snapshot_id", format!("sha256-{}", "a".repeat(64))),
        ("source_path", String::new()),
        ("source_path", "../genome.fna".to_owned()),
        ("source_path", "/etc/passwd".to_owned()),
        ("source_path", "https://example.org/genome.fna".to_owned()),
        ("source_path", "ncbi_dataset/data/../genome.fna".to_owned()),
        ("source_path", "missing.fna".to_owned()),
        ("source_path", "README.md".to_owned()),
        ("source_path", "🧬".repeat(5000)),
    ] {
        let mut invalid_request = request.clone();
        invalid_request[field] = json!(invalid);
        server.tool_error("prepared_fasta", invalid_request).await;
    }
    let prepared = server.successful("prepared_fasta", request).await;
    verify_prepared_paths(&store, &prepared);
    let opened = server
        .successful("dataset_open", json!({"path":"genes.csv"}))
        .await;
    assert_eq!(opened["row_count"], 5);
    server.finish().await;
}

async fn cli_prepared(store: &Path, request_file: &Path) -> std::process::Output {
    timeout(
        DEADLINE,
        Command::new(binary())
            .args(["prepared", "fasta"])
            .arg("--store-root")
            .arg(store)
            .arg("--request-file")
            .arg(request_file)
            .env("PATH", "")
            .env_remove("PYTHONPATH")
            .output(),
    )
    .await
    .unwrap()
    .unwrap()
}

#[tokio::test]
async fn one_shot_prepared_cli_reuses_native_snapshot_and_enforces_strict_bounded_json() {
    let fixture = Fixture::new();
    let store = prepare_refseq_fixture(&fixture);
    let request_file = fixture.root.join("request.json");
    fs::write(
        &request_file,
        serde_json::to_vec(&refseq_register_request()).unwrap(),
    )
    .unwrap();
    let registered =
        successful_cli_json(&cli_storage(&fixture, &store, "register", &request_file).await);
    fs::remove_dir_all(&fixture.data).unwrap();
    let request = prepared_request(registered["snapshot"]["snapshot_id"].as_str().unwrap());
    fs::write(&request_file, serde_json::to_vec(&request).unwrap()).unwrap();
    let prepared = successful_cli_json(&cli_prepared(&store, &request_file).await);
    assert_eq!(prepared["reused"], false);
    verify_prepared_paths(&store, &prepared);
    // A request exactly at the file-size limit is valid; bound applies to bytes.
    let mut at_limit = serde_json::to_vec(&request).unwrap();
    at_limit.resize(64 * 1024, b' ');
    fs::write(&request_file, at_limit).unwrap();
    let reused = successful_cli_json(&cli_prepared(&store, &request_file).await);
    assert_eq!(reused["reused"], true);
    assert_eq!(reused["recipe_id"], prepared["recipe_id"]);

    let mut extra = request.clone();
    extra["extra"] = json!(true);
    let mut path_escape = request.clone();
    path_escape["source_path"] = json!("../outside.fna");
    let mut overlong = request.clone();
    overlong["source_path"] = json!("🧬".repeat(5000));
    let mut missing = request.clone();
    missing.as_object_mut().unwrap().remove("snapshot_id");
    let duplicate = format!(
        "{{\"reference\":\"{PREPARED_REFERENCE}\",{}",
        &request.to_string()[1..]
    );
    for invalid in [
        "{}".to_owned(),
        "null".to_owned(),
        "[]".to_owned(),
        "{malformed JSON".to_owned(),
        missing.to_string(),
        duplicate,
        extra.to_string(),
        path_escape.to_string(),
        overlong.to_string(),
        format!("{request} {{}}"),
        " ".repeat(64 * 1024 + 1),
    ] {
        fs::write(&request_file, invalid).unwrap();
        let output = cli_prepared(&store, &request_file).await;
        assert!(!output.status.success());
        assert!(
            output.stdout.is_empty(),
            "errors must not look like JSON results"
        );
        assert!(!output.stderr.is_empty() && output.stderr.len() <= 2060);
    }
    fs::remove_file(&request_file).unwrap();
    let output = cli_prepared(&store, &request_file).await;
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    let output = cli_prepared(&store, &fixture.root).await;
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("request file must be a regular file"));
}

#[cfg(unix)]
#[tokio::test]
async fn prepared_cli_rejects_fifo_before_opening_it() {
    let fixture = Fixture::new();
    let store = prepare_refseq_fixture(&fixture);
    let request_file = fixture.root.join("prepared-request.fifo");
    assert!(std::process::Command::new("mkfifo")
        .arg(&request_file)
        .status()
        .unwrap()
        .success());
    let output = cli_prepared(&store, &request_file).await;
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("request file must be a regular file"));
}
