//! Router boundary tests. The installed wheel acceptance separately verifies a
//! real interpreter and RECORD pairing; the fake sibling here tests only argv,
//! stdio and process-status transport.
use std::{fs, process::Command};

fn binary() -> &'static str {
    env!("CARGO_BIN_EXE_biov")
}

#[test]
fn native_help_lists_one_public_command_and_explicit_mcp_routes() {
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
    let output = Command::new(binary()).arg("setup").output().unwrap();
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("unknown command"));
}

#[cfg(unix)]
#[test]
fn paired_bridge_preserves_literal_argv_stdio_and_exit_status() {
    use std::{io::Write, os::unix::fs::PermissionsExt, process::Stdio};
    let temp = tempfile::tempdir().unwrap();
    let executable = temp.path().join("biov");
    fs::copy(binary(), &executable).unwrap();
    let interpreter = temp.path().join("python");
    // Deliberately no Python behavior: this fixture verifies the exact explicit
    // bridge command and pipe inheritance, not interpreter/package pairing.
    fs::write(
        &interpreter,
        "#!/bin/sh\nprintf '%s\\n' \"$@\" >&2\ncat\nexit 37\n",
    )
    .unwrap();
    fs::set_permissions(&interpreter, fs::Permissions::from_mode(0o755)).unwrap();
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
        os::unix::{fs::PermissionsExt, process::ExitStatusExt},
        process::Stdio,
    };
    let temp = tempfile::tempdir().unwrap();
    let executable = temp.path().join("biov");
    fs::copy(binary(), &executable).unwrap();
    let interpreter = temp.path().join("python");
    fs::write(
        &interpreter,
        "#!/bin/sh\nprintf '%s\\n' \"$$\"\nexec /bin/sleep 30\n",
    )
    .unwrap();
    fs::set_permissions(&interpreter, fs::Permissions::from_mode(0o755)).unwrap();
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
    let directory = tempfile::tempdir().unwrap();
    let executable = directory.path().join("biov");
    fs::copy(binary(), &executable).unwrap();
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
