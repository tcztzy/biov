//! Thin local tool adapter. Native argv belongs to the selected upstream tool.
use biov_tools::ToolStore;
use std::{ffi::OsString, path::PathBuf, process::ExitCode};

pub fn run(args: impl IntoIterator<Item = OsString>) -> Result<ExitCode, String> {
    let mut args = args.into_iter();
    let action = args
        .next()
        .ok_or("missing tools action; expected setup, inspect or exec")?;
    if action == "--help" || action == "-h" {
        if args.next().is_some() {
            return Err("unexpected argument after help".into());
        }
        eprintln!("{}", HELP);
        return Ok(ExitCode::SUCCESS);
    }
    if action != "setup" && action != "inspect" && action != "exec" {
        return Err("unknown tools action; expected setup, inspect or exec".into());
    }
    let mut root = None;
    let mut pixi = None;
    let mut cwd = None;
    let mut no_install = false;
    let name = loop {
        let word = args
            .next()
            .ok_or("missing tool name; expected samtools or goatools")?;
        if word == "--no-install" && action == "exec" {
            if no_install {
                return Err("duplicate --no-install".into());
            }
            no_install = true;
            continue;
        }
        if word == "--environment-root" || word == "--pixi" || (word == "--cwd" && action == "exec")
        {
            let slot = if word == "--environment-root" {
                &mut root
            } else if word == "--pixi" {
                &mut pixi
            } else {
                &mut cwd
            };
            if slot.is_some() {
                return Err(format!("duplicate {}", word.to_string_lossy()));
            }
            let value = args.next().ok_or("option requires a path")?;
            if value.is_empty() || value.to_string_lossy().starts_with('-') {
                return Err("option requires a nonempty path".into());
            }
            *slot = Some(PathBuf::from(value));
            continue;
        }
        let name = word.to_str().ok_or("tool name must be UTF-8")?;
        if name.starts_with('-') {
            return Err(format!("unknown tools option {name}"));
        }
        break name.to_owned();
    };
    let mut native: Vec<_> = args.collect();
    if native.first().is_some_and(|arg| arg == "--") {
        native.remove(0);
    }
    if action != "exec" && !native.is_empty() {
        return Err("unexpected argument after tool name".into());
    }
    if let Some(cwd) = &cwd {
        if !cwd.is_dir() {
            return Err("execution cwd must be an existing directory".into());
        }
    }
    let store = ToolStore::from_environment(root, pixi)?;
    if action == "inspect" {
        let record = store.inspect(&name)?;
        println!(
            "{}",
            serde_json::to_string(&record).map_err(|e| e.to_string())?
        );
        return Ok(ExitCode::SUCCESS);
    }
    if action == "setup" {
        let record = store.setup(&name)?;
        eprintln!(
            "{}",
            serde_json::to_string(&record).map_err(|e| e.to_string())?
        );
        return Ok(ExitCode::SUCCESS);
    }
    if !no_install {
        store.setup(&name)?;
    }
    let status = store.execute(&name, &native, cwd.as_deref())?;
    #[cfg(unix)]
    {
        use std::os::unix::process::ExitStatusExt;
        if let Some(signal) = status.signal() {
            return Ok(ExitCode::from((128 + signal).min(255) as u8));
        }
    }
    Ok(ExitCode::from(status.code().unwrap_or(1) as u8))
}

pub const HELP: &str = "BioV native locked tools (initial Linux-64 slice)\n\nUsage:\n  biov-rs tools setup [--environment-root DIR] [--pixi FILE] samtools|goatools\n  biov-rs tools inspect [--environment-root DIR] [--pixi FILE] samtools|goatools\n  biov-rs tools exec [--no-install] [--environment-root DIR] [--pixi FILE] [--cwd DIR] samtools|goatools [ARGS]...\n\nOptions precede the tool name; all later arguments pass unchanged to the\nupstream entry point. BIOV_ENVIRONMENT_ROOT and BIOV_PIXI_BIN are supported.\nA matching existing Pixi 0.81.0 is required; no manager download, source switch,\nglobal PATH mutation, shell profile change or SSH execution occurs. Default exec\nperforms locked setup/repair when the recorded executable is unavailable. --no-install\nrequires that receipt and never provisions. Inspect reports recorded setup,\nnot independent package integrity or scientific validation. Native paths are literal;\nquoted ~ is not expanded. Pre-launch errors return 2; a tool may also return 2.\nUpdate/removal/cache\ncleanup, other tools and project manifests remain outside this migration slice.";
