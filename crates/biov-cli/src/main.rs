//! Thin native command-line adapter. MCP exclusively owns stdout while serving;
//! one-shot local commands write a single bounded JSON result there instead.
mod mcp;
mod model_cli;
mod prepared_cli;
mod python_bridge;
mod storage_cli;
mod tools_cli;

use std::{collections::BTreeMap, ffi::OsString, path::PathBuf, process::ExitCode};

use biov_data::DatasetStore;
use biov_prepared::PreparedStore;
use biov_storage::NativeStore;
use rmcp::{transport::stdio, ServiceExt};

const HELP: &str = "BioV: native lifecycle, storage and analysis with installed Python capabilities\n\nUsage:\n  biov install [options] <tool>\n  biov list [options]\n  biov uninstall [options] <tool>\n  biov tools inspect|exec ...\n  biov model download|inspect ...\n  biov storage register --store-root DIR --source-root DIR --request-file JSON\n  biov storage resolve --store-root DIR --request-file JSON\n  biov prepared fasta --store-root DIR --request-file JSON\n  biov mcp\n  biov mcp-native --data-root DIR --output-root DIR [--store-root DIR]\n  biov python <legacy command> ...\n  biov analyze|inspect-analysis|pixi|exec|run|update ...\n  biov [--config FILE] <Python-backed command> ...\n  biov --help\n  biov --version\n\ninstall, list, uninstall, tools, model, storage, prepared and mcp-native use Rust.\nThe existing mcp, analyze, inspect-analysis, pixi, exec, run and update\ncommands use the Python interpreter paired with this installed BioV package.\nLegacy setup is available as `biov python setup`.\nUse `<command> --help` for its options. Install BioV with `uv tool install biov`\nto make both runtimes available; native routing does not load Python; GOATOOLS installation invokes uv with the paired interpreter.\n\nBoth MCP routes use JSON-RPC over stdio. Native MCP requires explicit trusted\ndata/output roots and optionally a native snapshot store. Help, diagnostics and\nerrors go to stderr; MCP exclusively owns stdout while serving.";

#[derive(Debug, PartialEq, Eq)]
enum Command {
    Help,
    Version,
    Mcp {
        data_root: PathBuf,
        output_root: PathBuf,
        store_root: Option<PathBuf>,
    },
    StorageRegister {
        store_root: PathBuf,
        source_root: PathBuf,
        request_file: PathBuf,
    },
    StorageResolve {
        store_root: PathBuf,
        request_file: PathBuf,
    },
    PreparedFasta {
        store_root: PathBuf,
        request_file: PathBuf,
    },
}

fn parse_args(args: impl IntoIterator<Item = OsString>) -> Result<Command, String> {
    let mut args = args.into_iter();
    let command = args
        .next()
        .ok_or("missing command; expected mcp-native, storage or prepared")?;
    if command == "--help" || command == "-h" {
        if args.next().is_some() {
            return Err("unexpected argument after --help".into());
        }
        return Ok(Command::Help);
    }
    if command == "--version" || command == "-V" {
        if args.next().is_some() {
            return Err("unexpected argument after --version".into());
        }
        return Ok(Command::Version);
    }
    let (mode, allowed): (&str, &[&str]) = if command == "mcp-native" {
        (
            "mcp-native",
            &["--data-root", "--output-root", "--store-root"],
        )
    } else if command == "storage" {
        match args.next().as_deref() {
            Some(action) if action == "register" => (
                "register",
                &["--store-root", "--source-root", "--request-file"],
            ),
            Some(action) if action == "resolve" => ("resolve", &["--store-root", "--request-file"]),
            Some(action) if action == "--help" || action == "-h" => {
                if args.next().is_some() {
                    return Err("unexpected argument after --help".into());
                }
                return Ok(Command::Help);
            }
            Some(_) => return Err("unknown storage command; expected register or resolve".into()),
            None => return Err("missing storage command; expected register or resolve".into()),
        }
    } else if command == "prepared" {
        match args.next().as_deref() {
            Some(action) if action == "fasta" => ("fasta", &["--store-root", "--request-file"]),
            Some(action) if action == "--help" || action == "-h" => {
                if args.next().is_some() {
                    return Err("unexpected argument after --help".into());
                }
                return Ok(Command::Help);
            }
            Some(_) => return Err("unknown prepared command; expected fasta".into()),
            None => return Err("missing prepared command; expected fasta".into()),
        }
    } else {
        return Err("unknown command; expected mcp-native, storage or prepared".into());
    };
    let mut options = BTreeMap::new();
    while let Some(arg) = args.next() {
        if arg == "--help" || arg == "-h" {
            if args.next().is_some() {
                return Err("unexpected argument after --help".into());
            }
            return Ok(Command::Help);
        }
        let flag = allowed
            .iter()
            .copied()
            .find(|flag| arg == *flag)
            .ok_or_else(|| format!("unknown argument; expected {}", allowed.join(", ")))?;
        if options.contains_key(flag) {
            return Err(format!("duplicate {flag}"));
        }
        let kind = if flag == "--request-file" {
            "file"
        } else {
            "directory"
        };
        let value = args
            .next()
            .ok_or_else(|| format!("{flag} requires a {kind}"))?;
        if value.is_empty() || value.to_string_lossy().starts_with('-') {
            return Err(format!("{flag} requires a {kind}"));
        }
        options.insert(flag, PathBuf::from(value));
    }
    let mut required = |flag| {
        options
            .remove(flag)
            .ok_or_else(|| format!("missing required {flag}"))
    };
    match mode {
        "mcp-native" => Ok(Command::Mcp {
            data_root: required("--data-root")?,
            output_root: required("--output-root")?,
            store_root: options.remove("--store-root"),
        }),
        "register" => Ok(Command::StorageRegister {
            store_root: required("--store-root")?,
            source_root: required("--source-root")?,
            request_file: required("--request-file")?,
        }),
        "resolve" => Ok(Command::StorageResolve {
            store_root: required("--store-root")?,
            request_file: required("--request-file")?,
        }),
        "fasta" => Ok(Command::PreparedFasta {
            store_root: required("--store-root")?,
            request_file: required("--request-file")?,
        }),
        _ => unreachable!("parser modes are fixed above"),
    }
}

#[tokio::main(flavor = "current_thread")]
async fn main() -> ExitCode {
    let raw: Vec<_> = std::env::args_os().skip(1).collect();
    let first = raw.first().map(OsString::as_os_str);
    if first.is_some_and(|value| value == "tools") {
        return tool_result(tools_cli::run(raw.into_iter().skip(1)));
    }
    if first.is_some_and(|value| value == "model") {
        return tool_result(model_cli::run(raw.into_iter().skip(1)));
    }
    if first.is_some_and(|value| value == "install" || value == "list" || value == "uninstall") {
        return tool_result(tools_cli::run(raw));
    }
    if first.is_some_and(|value| value == "python") {
        return tool_result(python_bridge::run(raw.into_iter().skip(1)));
    }
    if first.is_some_and(python_bridge::is_python_route) {
        return tool_result(python_bridge::run(raw));
    }
    if raw.is_empty() {
        eprintln!("{HELP}");
        return ExitCode::from(2);
    }
    let command = match parse_args(std::env::args_os().skip(1)) {
        Ok(command) => command,
        Err(error) => {
            eprintln!("biov: {error}\n\n{HELP}");
            return ExitCode::from(2);
        }
    };
    let result = match command {
        Command::Help => {
            eprintln!("{HELP}");
            Ok(())
        }
        Command::Version => {
            eprintln!("biov {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        Command::StorageRegister {
            store_root,
            source_root,
            request_file,
        } => storage_cli::register(&store_root, &source_root, &request_file),
        Command::StorageResolve {
            store_root,
            request_file,
        } => storage_cli::resolve(&store_root, &request_file),
        Command::PreparedFasta {
            store_root,
            request_file,
        } => prepared_cli::fasta(&store_root, &request_file),
        Command::Mcp {
            data_root,
            output_root,
            store_root,
        } => serve_mcp(data_root, output_root, store_root).await,
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("biov: {}", mcp::bounded_text(&error, 2048));
            ExitCode::FAILURE
        }
    }
}

fn tool_result(result: Result<ExitCode, String>) -> ExitCode {
    match result {
        Ok(status) => status,
        Err(error) => {
            eprintln!("biov: {}", mcp::bounded_text(&error, 2048));
            ExitCode::from(2)
        }
    }
}

async fn serve_mcp(
    data_root: PathBuf,
    output_root: PathBuf,
    store_root: Option<PathBuf>,
) -> Result<(), String> {
    let store = DatasetStore::new(&data_root, &output_root).map_err(|error| error.to_string())?;
    let native_store = store_root
        .as_ref()
        .map(NativeStore::new)
        .transpose()
        .map_err(|error| error.to_string())?;
    let prepared_store = store_root
        .map(PreparedStore::new)
        .transpose()
        .map_err(|error| error.to_string())?;
    // The official SDK owns initialization, framing, routing, and EOF shutdown.
    // Never install a stdout logger: it would corrupt this transport.
    let service = mcp::DatasetServer::new(store, data_root, native_store, prepared_store)
        .serve(stdio())
        .await
        .map_err(|error| format!("MCP initialization failed: {error}"))?;
    service
        .waiting()
        .await
        .map_err(|error| format!("MCP service failed: {error}"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(words: &[&str]) -> Result<Command, String> {
        parse_args(words.iter().map(OsString::from))
    }

    #[test]
    fn parses_required_roots_in_either_order_and_optional_storage() {
        assert_eq!(
            args(&[
                "mcp-native",
                "--output-root",
                "results",
                "--data-root",
                "input"
            ])
            .unwrap(),
            Command::Mcp {
                data_root: "input".into(),
                output_root: "results".into(),
                store_root: None,
            }
        );
        assert_eq!(
            args(&[
                "mcp-native",
                "--store-root",
                "saved",
                "--data-root",
                "input",
                "--output-root",
                "results"
            ])
            .unwrap(),
            Command::Mcp {
                data_root: "input".into(),
                output_root: "results".into(),
                store_root: Some("saved".into()),
            }
        );
    }

    #[test]
    fn parses_storage_commands_without_a_source_for_resolution() {
        assert_eq!(
            args(&[
                "storage",
                "register",
                "--request-file",
                "request.json",
                "--source-root",
                "input",
                "--store-root",
                "saved"
            ])
            .unwrap(),
            Command::StorageRegister {
                store_root: "saved".into(),
                source_root: "input".into(),
                request_file: "request.json".into(),
            }
        );
        assert_eq!(
            args(&[
                "storage",
                "resolve",
                "--store-root",
                "saved",
                "--request-file",
                "request.json"
            ])
            .unwrap(),
            Command::StorageResolve {
                store_root: "saved".into(),
                request_file: "request.json".into(),
            }
        );
    }

    #[test]
    fn parses_prepared_fasta_without_an_original_source() {
        assert_eq!(
            args(&[
                "prepared",
                "fasta",
                "--request-file",
                "request.json",
                "--store-root",
                "saved"
            ])
            .unwrap(),
            Command::PreparedFasta {
                store_root: "saved".into(),
                request_file: "request.json".into()
            }
        );
    }

    #[test]
    fn rejects_unknown_missing_and_duplicate_arguments() {
        for invalid in [
            vec![],
            vec!["other"],
            vec!["mcp-native"],
            vec!["mcp-native", "--unknown"],
            vec!["mcp-native", "--data-root"],
            vec!["mcp-native", "--data-root", "--output-root", "results"],
            vec!["mcp-native", "--data-root", "-h"],
            vec![
                "mcp-native",
                "--data-root",
                "in",
                "--data-root",
                "in",
                "--output-root",
                "out",
            ],
            vec![
                "mcp-native",
                "--store-root",
                "saved",
                "--store-root",
                "saved",
            ],
            vec!["--help", "extra"],
            vec!["--version", "extra"],
            vec!["prepared"],
            vec!["prepared", "other"],
            vec!["prepared", "fasta"],
            vec!["prepared", "fasta", "--store-root", "saved"],
            vec!["prepared", "fasta", "--request-file", "request.json"],
            vec!["prepared", "fasta", "--source-root", "input"],
            vec![
                "prepared",
                "fasta",
                "--store-root",
                "saved",
                "--request-file",
                "r.json",
                "--request-file",
                "r.json",
            ],
            vec!["prepared", "--help", "extra"],
            vec!["storage"],
            vec!["storage", "other"],
            vec!["storage", "register"],
            vec!["storage", "resolve"],
            vec!["storage", "resolve", "--request-file"],
            vec!["storage", "resolve", "--request-file", ""],
            vec!["storage", "resolve", "--source-root", "input"],
            vec![
                "storage",
                "register",
                "--store-root",
                "saved",
                "--source-root",
                "input",
            ],
            vec![
                "storage",
                "register",
                "--store-root",
                "saved",
                "--request-file",
                "r.json",
            ],
            vec![
                "storage",
                "register",
                "--source-root",
                "input",
                "--request-file",
                "r.json",
            ],
            vec![
                "storage",
                "resolve",
                "--store-root",
                "saved",
                "--request-file",
                "r.json",
                "--request-file",
                "r.json",
            ],
            vec!["storage", "--help", "extra"],
        ] {
            assert!(args(&invalid).is_err(), "unexpectedly accepted {invalid:?}");
        }
    }

    #[test]
    fn exposes_help_and_version_without_roots() {
        for words in [
            &["--help"][..],
            &["mcp-native", "--help"],
            &["storage", "--help"],
            &["storage", "register", "--help"],
            &["storage", "resolve", "--help"],
            &["prepared", "--help"],
            &["prepared", "fasta", "--help"],
        ] {
            assert_eq!(args(words).unwrap(), Command::Help);
        }
        assert_eq!(args(&["--version"]).unwrap(), Command::Version);
    }
}
