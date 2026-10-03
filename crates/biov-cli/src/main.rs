//! Thin native command-line adapter. Stdout belongs exclusively to MCP.
mod mcp;

use std::{ffi::OsString, path::PathBuf, process::ExitCode};

use biov_data::DatasetStore;
use rmcp::{transport::stdio, ServiceExt};

const HELP: &str = "BioV native dataset MCP server\n\nUsage:\n  biov-rs mcp --data-root DIR --output-root DIR\n  biov-rs --help\n  biov-rs --version\n\nThe MCP server uses JSON-RPC over stdio. Input paths are restricted to\n--data-root; exports are restricted to --output-root. No Python runtime is used.\nHelp, diagnostics, and errors are written to stderr.";

#[derive(Debug, PartialEq, Eq)]
enum Command {
    Help,
    Version,
    Mcp {
        data_root: PathBuf,
        output_root: PathBuf,
    },
}

fn parse_args(args: impl IntoIterator<Item = OsString>) -> Result<Command, String> {
    let mut args = args.into_iter();
    match args.next().as_deref() {
        Some(command) if command == "--help" || command == "-h" => {
            if args.next().is_some() {
                return Err("unexpected argument after --help".into());
            }
            return Ok(Command::Help);
        }
        Some(command) if command == "--version" || command == "-V" => {
            if args.next().is_some() {
                return Err("unexpected argument after --version".into());
            }
            return Ok(Command::Version);
        }
        Some(command) if command == "mcp" => {}
        Some(_) => return Err("unknown command; expected mcp".into()),
        None => return Err("missing command; expected mcp".into()),
    }
    let mut data_root = None;
    let mut output_root = None;
    while let Some(arg) = args.next() {
        if arg == "--help" || arg == "-h" {
            if args.next().is_some() {
                return Err("unexpected argument after --help".into());
            }
            return Ok(Command::Help);
        }
        let (slot, flag) = if arg == "--data-root" {
            (&mut data_root, "--data-root")
        } else if arg == "--output-root" {
            (&mut output_root, "--output-root")
        } else {
            return Err("unknown argument; expected --data-root or --output-root".into());
        };
        if slot.is_some() {
            return Err(format!("duplicate {flag}"));
        }
        let value = args
            .next()
            .ok_or_else(|| format!("{flag} requires a directory"))?;
        if value.is_empty() || value.to_string_lossy().starts_with("--") {
            return Err(format!("{flag} requires a directory"));
        }
        *slot = Some(PathBuf::from(value));
    }
    Ok(Command::Mcp {
        data_root: data_root.ok_or("missing required --data-root DIR")?,
        output_root: output_root.ok_or("missing required --output-root DIR")?,
    })
}

#[tokio::main(flavor = "current_thread")]
async fn main() -> ExitCode {
    let command = match parse_args(std::env::args_os().skip(1)) {
        Ok(command) => command,
        Err(error) => {
            eprintln!("biov-rs: {error}\n\n{HELP}");
            return ExitCode::from(2);
        }
    };
    let (data_root, output_root) = match command {
        Command::Help => {
            eprintln!("{HELP}");
            return ExitCode::SUCCESS;
        }
        Command::Version => {
            eprintln!("biov-rs {}", env!("CARGO_PKG_VERSION"));
            return ExitCode::SUCCESS;
        }
        Command::Mcp {
            data_root,
            output_root,
        } => (data_root, output_root),
    };
    let store = match DatasetStore::new(&data_root, &output_root) {
        Ok(store) => store,
        Err(error) => {
            eprintln!("biov-rs: {error}");
            return ExitCode::FAILURE;
        }
    };
    // The official SDK owns initialization, framing, routing, and EOF shutdown.
    // Never install a stdout logger: it would corrupt this transport.
    let service = match mcp::DatasetServer::new(store).serve(stdio()).await {
        Ok(service) => service,
        Err(error) => {
            eprintln!("biov-rs: MCP initialization failed: {error}");
            return ExitCode::FAILURE;
        }
    };
    match service.waiting().await {
        Ok(_) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("biov-rs: MCP service failed: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(words: &[&str]) -> Result<Command, String> {
        parse_args(words.iter().map(OsString::from))
    }

    #[test]
    fn parses_required_roots_in_either_order() {
        assert_eq!(
            args(&["mcp", "--output-root", "results", "--data-root", "input"]).unwrap(),
            Command::Mcp {
                data_root: "input".into(),
                output_root: "results".into()
            }
        );
    }

    #[test]
    fn rejects_unknown_missing_and_duplicate_arguments() {
        for invalid in [
            vec![],
            vec!["other"],
            vec!["mcp"],
            vec!["mcp", "--unknown"],
            vec!["mcp", "--data-root"],
            vec!["mcp", "--data-root", "--output-root", "results"],
            vec![
                "mcp",
                "--data-root",
                "in",
                "--data-root",
                "in",
                "--output-root",
                "out",
            ],
            vec!["--help", "extra"],
        ] {
            assert!(args(&invalid).is_err(), "unexpectedly accepted {invalid:?}");
        }
    }

    #[test]
    fn exposes_help_and_version_without_roots() {
        assert_eq!(args(&["--help"]).unwrap(), Command::Help);
        assert_eq!(args(&["mcp", "--help"]).unwrap(), Command::Help);
        assert_eq!(args(&["--version"]).unwrap(), Command::Version);
    }
}
