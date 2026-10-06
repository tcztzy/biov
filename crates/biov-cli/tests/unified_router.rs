//! Router boundary tests. The installed wheel acceptance separately verifies a
//! real interpreter and RECORD pairing; the fake sibling here tests only argv,
//! stdio and process-status transport.
use std::process::Command;
#[cfg(unix)]
use std::{
    fs,
    sync::{Arc, Mutex, Weak},
};

#[cfg(unix)]
static ROUTER_FIXTURES: Mutex<Weak<[tempfile::TempDir; 3]>> = Mutex::new(Weak::new());

#[cfg(unix)]
fn prepare_router_fixtures() -> Arc<[tempfile::TempDir; 3]> {
    let mut shared = ROUTER_FIXTURES.lock().unwrap();
    if let Some(fixtures) = shared.upgrade() {
        return fixtures;
    }
    let fixtures = Arc::new({
        use std::os::unix::fs::PermissionsExt;
        let directories: [tempfile::TempDir; 3] =
            std::array::from_fn(|_| tempfile::tempdir().unwrap());
        for directory in &directories {
            fs::copy(env!("CARGO_BIN_EXE_biov"), directory.path().join("biov")).unwrap();
        }
        for (directory, script) in [
            (
                &directories[0],
                "#!/bin/sh\nprintf '%s\\n' \"$@\" >&2\ncat\nexit 37\n",
            ),
            (
                &directories[1],
                "#!/bin/sh\nprintf '%s\\n' \"$$\"\nexec /bin/sleep 30\n",
            ),
        ] {
            let interpreter = directory.path().join("python");
            fs::write(&interpreter, script).unwrap();
            fs::set_permissions(interpreter, fs::Permissions::from_mode(0o755)).unwrap();
        }
        directories
    });
    *shared = Arc::downgrade(&fixtures);
    fixtures
}

fn binary() -> &'static str {
    env!("CARGO_BIN_EXE_biov")
}

#[cfg(unix)]
#[test]
fn router_templates_are_complete_before_any_native_command_starts() {
    let directories = prepare_router_fixtures();
    let source = fs::read(binary()).unwrap();
    for directory in directories.iter() {
        assert_eq!(fs::read(directory.path().join("biov")).unwrap(), source);
    }
    assert!(directories[0].path().join("python").is_file());
    assert!(directories[1].path().join("python").is_file());
    assert!(!directories[2].path().join("python").exists());
}

#[test]
fn native_help_lists_one_public_command_and_explicit_mcp_routes() {
    #[cfg(unix)]
    let _fixtures = prepare_router_fixtures();
    let output = Command::new(binary()).arg("--help").output().unwrap();
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    let text = String::from_utf8(output.stderr).unwrap();
    for route in [
        "biov install",
        "biov list",
        "biov uninstall",
        "biov mcp",
        "biov mcp-native",
        "biov python",
    ] {
        assert!(text.contains(route), "missing {route}");
    }
    assert!(!text.contains("biov-rs"));
    assert!(!text.contains("biov tools setup"));
}

#[test]
fn legacy_mcp_does_not_become_native_or_use_path_python() {
    #[cfg(unix)]
    let _fixtures = prepare_router_fixtures();
    let output = Command::new(binary())
        .args(["mcp", "--data-root", "input", "--output-root", "output"])
        .env("PATH", "")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    let text = String::from_utf8(output.stderr).unwrap();
    assert!(text.contains("paired"));
    assert!(!text.contains("missing required --data-root"));
}

#[test]
fn removed_top_setup_is_rejected_without_starting_python() {
    #[cfg(unix)]
    let _fixtures = prepare_router_fixtures();
    let output = Command::new(binary()).arg("setup").output().unwrap();
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("unknown command"));
}

#[cfg(unix)]
#[test]
fn paired_bridge_preserves_literal_argv_stdio_and_exit_status() {
    use std::{io::Write, process::Stdio};
    let fixtures = prepare_router_fixtures();
    let temp = &fixtures[0];
    let executable = temp.path().join("biov");
    // Deliberately no Python behavior: this fixture verifies the exact explicit
    // bridge command and pipe inheritance, not interpreter/package pairing.
    let mut child = Command::new(&executable)
        .args([
            "python",
            "exec",
            "pypi:example",
            "",
            "O'Connor",
            "$(touch must-not-exist)",
            "--literal",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(b"inherited stdin\n")
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert_eq!(output.status.code(), Some(37));
    assert_eq!(output.stdout, b"inherited stdin\n");
    let argv = String::from_utf8(output.stderr).unwrap();
    assert_eq!(
        argv.lines().collect::<Vec<_>>(),
        [
            "-I",
            "-c",
            "import sys; sys.argv[0] = sys.argv.pop(1); from biov._bridge import main; main()",
            executable.to_str().unwrap(),
            "exec",
            "pypi:example",
            "",
            "O'Connor",
            "$(touch must-not-exist)",
            "--literal"
        ]
    );
    assert!(!temp.path().join("must-not-exist").exists());
}

#[cfg(unix)]
#[test]
fn python_bridge_replaces_process_and_receives_direct_signal() {
    use std::{
        io::{BufRead, BufReader},
        os::unix::process::ExitStatusExt,
        process::Stdio,
    };
    let fixtures = prepare_router_fixtures();
    let temp = &fixtures[1];
    let executable = temp.path().join("biov");
    let mut child = Command::new(&executable)
        .arg("mcp")
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut line = String::new();
    BufReader::new(child.stdout.take().unwrap())
        .read_line(&mut line)
        .unwrap();
    let interpreter_pid: u32 = line.trim().parse().unwrap();
    let replacement = interpreter_pid == child.id();
    // Always clean up the interpreter even if a regression restores a wrapper.
    assert!(Command::new("/bin/kill")
        .args(["-TERM", &interpreter_pid.to_string()])
        .status()
        .unwrap()
        .success());
    let status = child.wait().unwrap();
    assert!(
        replacement,
        "Python bridge must replace the original process, not orphan a child"
    );
    assert_eq!(status.signal(), Some(15));
}

#[test]
fn lifecycle_rejects_removed_custom_publication_options_without_side_effects() {
    #[cfg(unix)]
    let _fixtures = prepare_router_fixtures();
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path().join("environments");
    for args in [
        vec!["list", "--json"],
        vec!["install", "--bin-dir", "/unrequested-bin", "samtools"],
        vec!["installed-run", "samtools"],
        vec!["install", "samtools", "--uv", "uv"],
    ] {
        let output = Command::new(binary())
            .args(args)
            .env("BIOV_ENVIRONMENT_ROOT", &root)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(!root.exists());
    }
}

#[test]
fn empty_global_list_requires_no_managers_or_directory_creation() {
    #[cfg(unix)]
    let _fixtures = prepare_router_fixtures();
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path().join("environments");
    let output = Command::new(binary())
        .args(["list", "--environment-root", root.to_str().unwrap()])
        .env("PATH", "")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stdout).contains("No"));
    assert!(!root.exists());
}

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
#[test]
fn goatools_install_requires_paired_python_before_manager_or_root_mutation() {
    let fixtures = prepare_router_fixtures();
    let directory = &fixtures[2];
    let executable = directory.path().join("biov");
    let root = directory.path().join("environments");
    let output = Command::new(executable)
        .args([
            "install",
            "--environment-root",
            root.to_str().unwrap(),
            "--uv",
            "/must-not-run-uv",
            "goatools",
        ])
        .env("PATH", "")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("paired"));
    assert!(!root.exists());
}
