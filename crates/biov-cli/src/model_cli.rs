//! Thin selected-file model download and offline verification adapter.
use biov_tools::model_resources::{self, ModelDownload};
use std::{ffi::OsString, path::PathBuf, process::ExitCode};

pub fn run(arguments: impl IntoIterator<Item = OsString>) -> Result<ExitCode, String> {
    let mut arguments = arguments.into_iter();
    let action = arguments
        .next()
        .ok_or("missing model action; expected download or inspect")?;
    if action == "--help" || action == "-h" {
        if arguments.next().is_some() {
            return Err("unexpected argument after help".into());
        }
        eprintln!("{HELP}");
        return Ok(ExitCode::SUCCESS);
    }
    if action == "inspect" {
        let directory = arguments
            .next()
            .ok_or("model inspect requires a resource directory")?;
        if directory == "--help" || directory == "-h" {
            if arguments.next().is_some() {
                return Err("unexpected argument after help".into());
            }
            eprintln!("{HELP}");
            return Ok(ExitCode::SUCCESS);
        }
        if directory.is_empty()
            || directory.to_string_lossy().starts_with('-')
            || arguments.next().is_some()
        {
            return Err("model inspect requires exactly one directory".into());
        }
        let result = model_resources::inspect(&PathBuf::from(directory))?;
        println!(
            "{}",
            serde_json::to_string(&result).map_err(|e| e.to_string())?
        );
        return Ok(ExitCode::SUCCESS);
    }
    if action != "download" {
        return Err("unknown model action; expected download or inspect".into());
    }
    let mut revision = None;
    let mut local_dir = None;
    let mut hf = None;
    let mut uv = None;
    let mut no_install = false;
    let repository = loop {
        let word = arguments
            .next()
            .ok_or("model download requires a repository and explicit filenames")?;
        if word == "--help" || word == "-h" {
            if arguments.next().is_some() {
                return Err("unexpected argument after help".into());
            }
            eprintln!("{HELP}");
            return Ok(ExitCode::SUCCESS);
        }
        if word == "--no-install" {
            if no_install {
                return Err("duplicate --no-install".into());
            }
            no_install = true;
            continue;
        }
        if word == "--revision" || word == "--local-dir" || word == "--hf" || word == "--uv" {
            let slot = if word == "--revision" {
                &mut revision
            } else if word == "--local-dir" {
                &mut local_dir
            } else if word == "--hf" {
                &mut hf
            } else {
                &mut uv
            };
            if slot.is_some() {
                return Err(format!("duplicate {}", word.to_string_lossy()));
            }
            let value = arguments.next().ok_or("option requires a value")?;
            if value.is_empty() || value.to_string_lossy().starts_with('-') {
                return Err("option requires a nonempty value".into());
            }
            *slot = Some(value);
            continue;
        }
        if word.to_string_lossy().starts_with('-') {
            return Err("unknown model download option".into());
        }
        break word
            .into_string()
            .map_err(|_| "model repository must be UTF-8")?;
    };
    let files = arguments
        .map(|value| {
            value
                .into_string()
                .map_err(|_| "model filenames must be UTF-8".to_owned())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let request = ModelDownload {
        repository,
        revision: revision
            .ok_or("missing --revision; provide a full immutable Git commit")?
            .into_string()
            .map_err(|_| "revision must be UTF-8")?,
        files,
        local_dir: PathBuf::from(local_dir.ok_or("missing --local-dir")?),
        hf: hf.map(PathBuf::from),
        uv: uv.map(PathBuf::from),
        no_install,
    };
    let result = model_resources::download(&request)?;
    println!(
        "{}",
        serde_json::to_string(&result).map_err(|e| e.to_string())?
    );
    Ok(ExitCode::SUCCESS)
}

pub const HELP: &str = "BioV model resources: selected native files, official hf downloads\n\nUsage:\n  biov model download --revision COMMIT --local-dir DIR [--hf FILE] [--uv FILE] [--no-install] REPO FILE...\n  biov model inspect DIR\n\nOptions precede the repository. COMMIT must be an immutable full Git commit;\nfilenames must be explicit normalized relative files, without globs or folders.\nA new directory is published after official hf succeeds and every selected file\nis verified. Existing matching complete resources are reused offline; conflicting\nor damaged directories fail unchanged. Inspect hashes complete selected files\nwithout hf, uv, Python, network or model execution. This is selected-file storage,\nnot a claim that a complete or usable model was downloaded.\n\nA compatible installed hf on PATH is reused; --hf/BIOV_HF_BIN selects it explicitly.\nWhen none is compatible, uv tool run supplies huggingface-hub==2.1.1 in its native\non-demand environment, without downloading Python. --uv/BIOV_UV_BIN chooses uv.\n--no-install disables fallback provisioning, not resource network access.\nBioV adds no token arguments and records no credentials; official hf owns existing\nauthentication, network behavior and cache metadata. Custom HF_ENDPOINT mirrors\nare outside this initial source contract. Failed unpublished downloads are retained\nwith their location in diagnostics. No remove/update/cache-cleanup command is\nimplemented. Linux x86_64 download publication is the initial supported platform.\nNative payloads plus JSON/checksum/README companions remain readable without BioV.\nSuccessful commands emit one JSON result on stdout; upstream progress and errors\nuse stderr. Errors return 2; the original upstream status is reported in diagnostics.";
