//! Portable descriptions complement (and never weaken) the strict reopen pair.
use biov_data::{DatasetStore, ExportRequest, OpenRequest, QueryRequest, ReopenRequest};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{fs, path::Path};

fn setup(csv: &str, metadata: Value) -> (tempfile::TempDir, DatasetStore, Value) {
    let dir = tempfile::tempdir().unwrap();
    fs::create_dir(dir.path().join("input")).unwrap();
    fs::create_dir(dir.path().join("output")).unwrap();
    fs::write(dir.path().join("input/table.csv"), csv).unwrap();
    let mut store =
        DatasetStore::new(&dir.path().join("input"), &dir.path().join("output")).unwrap();
    let opened = store
        .open(
            serde_json::from_value::<OpenRequest>(json!({
                "path":"table.csv", "preview_rows":0, "metadata":metadata,
                "schema":{"count":"int64","ratio":"float64","passed":"boolean"}
            }))
            .unwrap(),
        )
        .unwrap();
    (dir, store, opened)
}

fn export(store: &mut DatasetStore, opened: &Value) -> Value {
    store
        .export(ExportRequest {
            dataset_id: opened["dataset_id"].as_str().unwrap().into(),
        })
        .unwrap()
}

fn load_manifest(exported: &Value) -> Value {
    serde_json::from_slice(&fs::read(exported["manifest_path"].as_str().unwrap()).unwrap()).unwrap()
}

#[test]
fn complete_bundle_has_predictable_relative_files_and_content_identities() {
    let csv = "id,count,ratio,passed\n001,9007199254740993,1.5,true\n002,,,\n";
    let (dir, mut store, opened) = setup(csv, json!({}));
    let exported = export(&mut store, &opened);
    let manifest = load_manifest(&exported);
    let artifact = exported["artifact_id"].as_str().unwrap();
    assert_eq!(manifest["manifest_version"], 1);
    assert_eq!(manifest["format"], "arrow_ipc_file");
    assert_eq!(manifest["artifact_id"], artifact);
    assert_eq!(manifest["files"]["arrow"], format!("{artifact}.arrow"));
    assert_eq!(manifest["files"]["record"], format!("{artifact}.json"));
    assert_eq!(manifest["files"]["readme"], format!("{artifact}.README.md"));
    assert_eq!(
        Path::new(exported["manifest_path"].as_str().unwrap())
            .file_name()
            .unwrap(),
        format!("{artifact}.manifest.json").as_str()
    );
    for filename in manifest["files"].as_object().unwrap().values() {
        let relative = Path::new(filename.as_str().unwrap());
        assert!(!relative.is_absolute());
        assert_eq!(relative.components().count(), 1);
        assert!(dir.path().join("output").join(relative).is_file());
    }
    assert_eq!(fs::read_dir(dir.path().join("output")).unwrap().count(), 4);
    for (path, identity) in [
        ("execution_host_path", "content"),
        ("record_path", "record"),
    ] {
        let bytes = fs::read(exported[path].as_str().unwrap()).unwrap();
        assert_eq!(manifest[identity]["bytes"], bytes.len());
        assert_eq!(
            manifest[identity]["sha256"],
            format!("{:x}", Sha256::digest(&bytes))
        );
    }
    assert_eq!(manifest["content"]["row_count"], 2);
    assert_eq!(manifest["record"]["record_version"], 2);
    assert_eq!(manifest["source"]["bytes"], csv.len());
    assert_eq!(
        manifest["source"]["sha256"],
        format!("{:x}", Sha256::digest(csv))
    );
    assert_eq!(manifest["source"]["historical_path"], "table.csv");
    assert_eq!(manifest["source"]["required_for_reading"], false);
    assert_eq!(
        manifest["source"]["path_interpretation"],
        "relative_to_original_data_root_not_bundle"
    );
    // Large companions stay in files rather than expanding the MCP response.
    assert!(exported.get("manifest").is_none());
    assert!(exported.get("readme").is_none());
    assert!(
        fs::metadata(exported["manifest_path"].as_str().unwrap())
            .unwrap()
            .len()
            <= 128 * 1024
    );
    assert!(
        fs::metadata(exported["readme_path"].as_str().unwrap())
            .unwrap()
            .len()
            <= 16 * 1024
    );
}

#[test]
fn dictionary_is_typed_ordered_and_preserves_unknown_scientific_semantics() {
    let (_dir, mut store, opened) = setup(
        "id,count,ratio,passed\n001,10,1.5,true\n002,,,\n",
        json!({"species":"Homo sapiens","reference":"GRCh38","units":"reads","coordinates":"1-based closed"}),
    );
    let manifest = load_manifest(&export(&mut store, &opened));
    for (column, (name, logical_type, null_count)) in
        manifest["columns"].as_array().unwrap().iter().zip([
            ("id", "string", 0),
            ("count", "int64", 1),
            ("ratio", "float64", 1),
            ("passed", "boolean", 1),
        ])
    {
        assert_eq!(column["name"], name);
        assert_eq!(column["logical_type"], logical_type);
        assert!(column["polars_dtype"].is_string());
        assert_eq!(column["nullable"], true);
        assert_eq!(column["null_count"], null_count);
        for unknown in ["description", "units", "coordinates"] {
            assert!(column[unknown].is_null());
        }
    }
    assert_eq!(manifest["scientific_metadata"]["reference"], "GRCh38");
    assert_eq!(manifest["scientific_metadata"]["units"], "reads");
    assert_eq!(manifest["biological_identifier"], Value::Null);
    assert_eq!(manifest["versions"]["reference_version"], Value::Null);
    assert_eq!(manifest["versions"]["provider_release"], Value::Null);
    assert_eq!(manifest["trust"]["authenticity"], "not_established");
}

#[test]
fn biological_accession_and_provider_versions_have_distinct_meanings() {
    for (declared, namespace, accession, base, version, kind) in [
        (
            "refseq.gcf://GCF_000001405.040",
            "refseq.gcf",
            "GCF_000001405.040",
            "GCF_000001405",
            json!("040"),
            json!("assembly_revision"),
        ),
        (
            "GCF_000001405",
            "refseq.gcf",
            "GCF_000001405",
            "GCF_000001405",
            Value::Null,
            json!("assembly_revision"),
        ),
        (
            "uniprot:P05067",
            "uniprot",
            "P05067",
            "P05067",
            Value::Null,
            Value::Null,
        ),
    ] {
        let (_dir, mut store, opened) = setup(
            "id,count,ratio,passed\n001,10,1.5,true\n",
            json!({"identifier":declared}),
        );
        let manifest = load_manifest(&export(&mut store, &opened));
        let identifier = &manifest["biological_identifier"];
        assert_eq!(identifier["namespace"], namespace);
        assert_eq!(identifier["accession"], accession);
        assert_eq!(identifier["base_accession"], base);
        assert_eq!(identifier["accession_version"], version);
        assert_eq!(identifier["accession_version_kind"], kind);
        for key in ["entry_version", "sequence_version", "provider_release"] {
            assert_eq!(identifier[key], Value::Null);
        }
        assert_eq!(
            identifier["validation"],
            "syntax_only_not_provider_verified"
        );
        assert_eq!(manifest["scientific_metadata"]["identifier"], declared);
        assert!(manifest["version_semantics"]["record_version"].is_string());
    }
}

#[test]
fn moved_bundle_keeps_lineage_without_source_files_or_old_session() {
    let (dir, mut store, opened) = setup(
        "id,count,ratio,passed\n001,10,1.5,true\n002,30,2.5,false\n003,20,,\n",
        json!({}),
    );
    let derived = store
        .query(
            serde_json::from_value::<QueryRequest>(json!({
                "dataset_id":opened["dataset_id"],
                "filter":{"column":"count","op":"ge","value":20},
                "sort":{"column":"count","descending":true},
                "select":["id","count"], "preview_rows":0
            }))
            .unwrap(),
        )
        .unwrap();
    let exported = export(&mut store, &derived);
    let manifest = load_manifest(&exported);
    assert_eq!(
        manifest["lineage"]["operations"].as_array().unwrap().len(),
        1
    );
    assert_eq!(manifest["lineage"]["historical_references_only"], true);
    assert_eq!(manifest["lineage"]["reopen_verification"], Value::Null);
    assert_eq!(manifest["columns"].as_array().unwrap().len(), 2);
    let moved = tempfile::tempdir().unwrap();
    for entry in fs::read_dir(dir.path().join("output")).unwrap() {
        let entry = entry.unwrap();
        fs::copy(entry.path(), moved.path().join(entry.file_name())).unwrap();
    }
    drop(store);
    drop(dir);
    let mut restarted = DatasetStore::new(moved.path(), moved.path()).unwrap();
    let reopened = restarted
        .reopen(ReopenRequest {
            record_path: manifest["files"]["record"].as_str().unwrap().into(),
            preview_rows: 5,
        })
        .unwrap();
    assert_eq!(
        reopened["preview"]["rows"],
        json!([["002", 30], ["003", 20]])
    );
    let next = export(&mut restarted, &reopened);
    let next_manifest = load_manifest(&next);
    assert_eq!(
        next_manifest["lineage"]["operations"],
        manifest["lineage"]["operations"]
    );
    assert_eq!(
        next_manifest["lineage"]["reopen_verification"],
        reopened["reopen_verification"]
    );
    assert_eq!(next_manifest["source"], manifest["source"]);
}

#[test]
fn caller_metadata_is_escaped_json_data_never_interpolated_into_readme_code() {
    let dir = tempfile::tempdir().unwrap();
    let name = "\"\\```python [column](https://example.invalid) 基因";
    let declaration = "\"\\```python print('UNTRUSTED_METADATA')";
    let mut csv = csv::Writer::from_writer(Vec::new());
    csv.write_record([name]).unwrap();
    csv.write_record(["0001"]).unwrap();
    fs::write(dir.path().join("table.csv"), csv.into_inner().unwrap()).unwrap();
    let mut store = DatasetStore::new(dir.path(), dir.path()).unwrap();
    let opened = store
        .open(
            serde_json::from_value(json!({
                "path":"table.csv", "metadata":{"species":declaration}
            }))
            .unwrap(),
        )
        .unwrap();
    let exported = export(&mut store, &opened);
    let manifest = load_manifest(&exported);
    assert_eq!(manifest["columns"][0]["name"], name);
    assert_eq!(manifest["scientific_metadata"]["species"], declaration);
    let readme = fs::read_to_string(exported["readme_path"].as_str().unwrap()).unwrap();
    assert!(!readme.contains(name));
    assert!(!readme.contains(declaration));
    assert!(!readme.contains("UNTRUSTED_METADATA"));
    assert_eq!(readme.matches("```python").count(), 1);
    assert_eq!(readme.matches("```").count(), 2);
    assert!(readme.contains("ipc.open_file"));
    assert!(readme.contains("pc.is_valid"));
    assert!(readme.contains("pc.mean"));
}

#[test]
fn descriptive_companions_are_not_a_new_requirement_for_strict_reopen() {
    let (dir, mut store, opened) = setup("id,count,ratio,passed\n001,10,1.5,true\n", json!({}));
    let exported = export(&mut store, &opened);
    let manifest = load_manifest(&exported);
    fs::remove_file(exported["manifest_path"].as_str().unwrap()).unwrap();
    fs::write(
        exported["readme_path"].as_str().unwrap(),
        "untrusted changes",
    )
    .unwrap();
    let output = dir.path().join("output");
    let reopened = DatasetStore::new(&output, &output)
        .unwrap()
        .reopen(ReopenRequest {
            record_path: manifest["files"]["record"].as_str().unwrap().into(),
            preview_rows: 1,
        })
        .unwrap();
    assert_eq!(reopened["preview"]["rows"], json!([["001", 10, 1.5, true]]));
}
