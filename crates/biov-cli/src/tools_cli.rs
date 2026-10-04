//! Thin local tool adapter. Native argv belongs to the selected upstream tool.
use biov_tools::ToolStore;
use std::{ffi::OsString, path::PathBuf, process::ExitCode};

pub fn run(args: impl IntoIterator<Item = OsString>) -> Result<ExitCode, String> {
    let mut args = args.into_iter();
    let action = args
        .next()
        .ok_or("missing tools action; expected inspect or exec")?;
    if action == "--help" || action == "-h" {
        if args.next().is_some() {
            return Err("unexpected argument after help".into());
        }
        eprintln!("{}", HELP);
        return Ok(ExitCode::SUCCESS);
    }
    if action == "install" || action == "list" || action == "uninstall" {
        return installed_cli(&action, args);
    }
    if action != "inspect" && action != "exec" {
        return Err("unknown tools action; expected inspect or exec".into());
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
    if !no_install {
        store.setup(&name)?;
    }
    let status = store.execute_with_interrupt_forwarding(&name, &native, cwd.as_deref())?;
    #[cfg(unix)]
    {
        use std::os::unix::process::ExitStatusExt;
        if let Some(signal) = status.signal() {
            return Ok(ExitCode::from((128 + signal).min(255) as u8));
        }
    }
    Ok(ExitCode::from(status.code().unwrap_or(1) as u8))
}

pub const HELP: &str = "BioV native locked tools (initial Linux-64 slice)\n\nUsage:\n  biov tools inspect [--environment-root DIR] [--pixi FILE] samtools|goatools\n  biov tools exec [--no-install] [--environment-root DIR] [--pixi FILE] [--cwd DIR] samtools|goatools [ARGS]...\n\nOptions precede the tool name; all later arguments pass unchanged to the\nupstream entry point. BIOV_ENVIRONMENT_ROOT and BIOV_PIXI_BIN are supported.\nA matching existing Pixi 0.81.0 is required; no manager download, source switch,\nglobal PATH mutation, shell profile change or SSH execution occurs. Default exec\nperforms locked setup/repair when the recorded executable is unavailable. --no-install\nrequires that receipt and never provisions. Inspect reports recorded setup,\nnot independent package integrity or scientific validation. Native paths are literal;\nquoted ~ is not expanded. Pre-launch errors return 2; a tool may also return 2.\nSoftware upgrades/cache\ncleanup, other tools and project manifests remain outside this migration slice.";

fn installed_cli(
    action: &OsString,
    args: impl IntoIterator<Item = OsString>,
) -> Result<ExitCode, String> {
    let mut args = args.into_iter();
    let mut root = None;
    let mut pixi = None;
    let mut uv = None;
    let mut name = None;
    while let Some(word) = args.next() {
        if word == "--help" || word == "-h" {
            eprintln!("{INSTALL_HELP}");
            return Ok(ExitCode::SUCCESS);
        }
        if word == "--environment-root" || word == "--pixi" || word == "--uv" {
            if name.is_some() {
                return Err("options must precede the tool name".into());
            }
            let slot = if word == "--environment-root" {
                &mut root
            } else if word == "--pixi" {
                &mut pixi
            } else {
                &mut uv
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
        if word.to_string_lossy().starts_with('-') || action == "list" || name.is_some() {
            return Err(format!("unexpected argument {}", word.to_string_lossy()));
        }
        name = Some(word.to_str().ok_or("tool name must be UTF-8")?.to_owned());
    }
    let store = ToolStore::from_environment(root, pixi)?;
    if action == "list" {
        print!("{}", store.list_global(uv.as_deref())?);
    } else {
        let name = name.ok_or("missing tool name; expected samtools or goatools")?;
        if action == "uninstall" {
            store.uninstall_global(&name, uv.as_deref())?;
        } else {
            let python = if name == "goatools" {
                Some(crate::python_bridge::installed_interpreter()?)
            } else {
                None
            };
            store.install_global(&name, python.as_deref(), uv.as_deref())?;
            for directory in store.bin_paths() {
                if !std::env::split_paths(&std::env::var_os("PATH").unwrap_or_default()).any(|p| {
                    p == directory
                        || std::fs::canonicalize(p).ok().as_deref() == Some(directory.as_path())
                }) {
                    let p = directory
                        .to_str()
                        .ok_or("bin path must be UTF-8")?
                        .replace('\'', "'\\''");
                    eprintln!(
                        "Tool bin is not on PATH. For sh/bash/zsh, run: export PATH='{}':\"$PATH\"\nFor fish, run: fish_add_path '{}'\nAdd the matching command to your shell configuration for future sessions. BioV did not edit any shell profile.",
                        p, p
                    );
                }
            }
        }
    }
    Ok(ExitCode::SUCCESS)
}

pub const INSTALL_HELP: &str = "Installed user software (Linux x86_64)\n\nUsage:\n  biov install [--environment-root DIR] [--pixi FILE] [--uv FILE] samtools|goatools\n  biov list [--environment-root DIR] [--pixi FILE] [--uv FILE]\n  biov uninstall [--environment-root DIR] [--pixi FILE] [--uv FILE] samtools|goatools\n\nPixi global manages Samtools; uv tool manages GOATOOLS with the Python interpreter\npaired with this BioV installation. BIOV_UV_BIN selects uv. Managers own isolated\nenvironments and command bins under the BioV environment root. Add the printed\nbin paths to PATH; no shell profile or user global settings are changed.\nList shows the upstream inventory only in these dedicated roots. No JSON\ninterface is claimed for uv's human-readable list. Uninstall delegates to the\nbackend and removes the selected user-tool environment and commands; scientific\ndata and separate locked-workflow environments are not removed.\nGlobal installs pin primary package versions but do not consume BioV's full\nscientific lock. Use tools exec for the separate bundled locked workflow.";
