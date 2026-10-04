//! Explicit compatibility boundary for Python-backed capabilities. Native routes
//! never load Python; installed legacy routes replace this process with the
//! interpreter beside the real installed executable, never an ambient PATH one.
use std::{ffi::OsStr, ffi::OsString, path::Path, process::Command, process::ExitCode};

// Set the legacy argv before importing biov: its configuration intentionally
// preselects --config during package import, before Typer runs its callback.
const BOOTSTRAP: &str =
    "import sys; sys.argv[0] = sys.argv.pop(1); from biov._bridge import main; main()";

pub fn is_python_route(value: &OsStr) -> bool {
    value
        .to_str()
        .is_some_and(|word| word.starts_with("--config="))
        || matches!(
            value.to_str(),
            Some(
                "mcp"
                    | "analyze"
                    | "inspect-analysis"
                    | "pixi"
                    | "exec"
                    | "run"
                    | "update"
                    | "--config"
            )
        )
}

fn paired_interpreter(binary: &Path) -> Result<std::path::PathBuf, String> {
    let directory = binary
        .parent()
        .ok_or("installed executable has no parent")?;
    let interpreter = directory.join(if cfg!(windows) {
        "python.exe"
    } else {
        "python"
    });
    if !interpreter.is_file() {
        return Err(format!(
            "this Python-backed command needs the interpreter paired with the installed BioV package; none exists at {}. Install BioV with `uv tool install biov`. No Python executable from PATH is used",
            interpreter.display()
        ));
    }
    Ok(interpreter)
}

pub fn installed_interpreter() -> Result<std::path::PathBuf, String> {
    let binary = std::env::current_exe()
        .and_then(std::fs::canonicalize)
        .map_err(|error| format!("cannot locate installed BioV executable: {error}"))?;
    paired_interpreter(&binary)
}

pub fn run(arguments: impl IntoIterator<Item = OsString>) -> Result<ExitCode, String> {
    // current_exe returns the real binary on supported systems, including uv's
    // public tool symlink. Do not canonicalize the interpreter itself: its venv
    // symlink path is how Python selects the paired environment.
    let binary = std::env::current_exe()
        .and_then(std::fs::canonicalize)
        .map_err(|error| format!("cannot locate installed BioV executable: {error}"))?;
    let interpreter = paired_interpreter(&binary)?;
    let mut command = Command::new(interpreter);
    command
        .arg("-I")
        .arg("-c")
        .arg(BOOTSTRAP)
        .arg(&binary)
        .args(arguments);
    // -I excludes cwd, PYTHONPATH and user site packages. The Python bridge also
    // checks the distribution RECORD, so a stale checkout/package cannot supply
    // this installed binary's implementation. All stdio and argv are inherited.
    #[cfg(unix)]
    {
        use std::os::unix::process::CommandExt;
        let error = command.exec();
        Err(format!("cannot start paired Python interpreter: {error}"))
    }
    #[cfg(not(unix))]
    {
        let status = command
            .status()
            .map_err(|error| format!("cannot start paired Python interpreter: {error}"))?;
        // Windows process exit codes are wider than ExitCode::from(u8).
        // Preserve the child status rather than truncating it to one byte.
        std::process::exit(status.code().unwrap_or(1));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_routes_are_explicit_and_native_routes_are_excluded() {
        for route in [
            "mcp",
            "analyze",
            "pixi",
            "exec",
            "run",
            "update",
            "--config",
            "--config=config.toml",
        ] {
            assert!(is_python_route(OsStr::new(route)), "{route}");
        }
        for route in [
            "mcp-native",
            "storage",
            "prepared",
            "tools",
            "install",
            "list",
            "uninstall",
            "setup",
            "--help",
            "--version",
            "unknown",
        ] {
            assert!(!is_python_route(OsStr::new(route)), "{route}");
        }
    }

    #[test]
    fn no_ambient_interpreter_fallback() {
        let directory = tempfile::tempdir().unwrap();
        assert!(paired_interpreter(&directory.path().join("biov")).is_err());
    }
}
