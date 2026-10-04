#![cfg(all(target_os = "linux", target_arch = "x86_64"))]
use std::{fs, process::Command};
#[path = "../../biov-tools/tests/support/mod.rs"]
mod support;

#[test]
fn native_cli_setup_inspection_literal_argv_and_status() {
    let dir = tempfile::tempdir().unwrap();
    let pixi = support::manager(dir.path());
    let root = dir.path().join("environments");
    let run = |words: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_biov-rs"))
            .args(words)
            .env("BIOV_ENVIRONMENT_ROOT", &root)
            .env("BIOV_PIXI_BIN", &pixi)
            .output()
            .unwrap()
    };
    let before = run(&["tools", "inspect", "goatools"]);
    assert!(before.status.success());
    assert_eq!(
        serde_json::from_slice::<serde_json::Value>(&before.stdout).unwrap()["status"],
        "unavailable"
    );
    assert!(!root.exists());
    let unavailable_cwd = dir.path().join("unavailable");
    let bad_cwd = run(&[
        "tools",
        "exec",
        "--cwd",
        unavailable_cwd.to_str().unwrap(),
        "goatools",
    ]);
    assert_eq!(bad_cwd.status.code(), Some(2));
    assert!(!root.exists());
    let missing = run(&["tools", "exec", "--no-install", "goatools", "--help"]);
    assert_eq!(missing.status.code(), Some(2));
    assert!(!root.exists());
    let first = run(&[
        "tools",
        "exec",
        "--cwd",
        dir.path().to_str().unwrap(),
        "goatools",
        "--",
        "",
        "O'Connor",
        "$(touch must-not-exist)",
        "--environment-root",
        "native-option",
    ]);
    assert_eq!(
        first.status.code(),
        Some(37),
        "{}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert_eq!(first.stdout, b"native stdout\n");
    assert!(String::from_utf8(first.stderr)
        .unwrap()
        .contains("native stderr"));
    let native: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.path().join("native.json")).unwrap()).unwrap();
    assert_eq!(native["argv"][1], "");
    assert_eq!(native["argv"][2], "O'Connor");
    assert_eq!(native["argv"][3], "$(touch must-not-exist)");
    assert_eq!(native["argv"][4], "--environment-root");
    let second = run(&["tools", "exec", "--no-install", "goatools", "--help"]);
    assert_eq!(second.status.code(), Some(37));
    let log = fs::read_to_string(dir.path().join("argv.jsonl")).unwrap();
    assert_eq!(
        log.lines()
            .filter(|line| line.starts_with("[\"install\""))
            .count(),
        1
    );
    assert!(!dir.path().join("must-not-exist").exists());
}

#[test]
fn parser_errors_and_help_have_no_execution_side_effects() {
    // Synchronize fixture materialization before either test starts a child.
    support::prepare();
    for args in [
        vec!["tools"],
        vec!["tools", "unknown"],
        vec!["tools", "exec", "--pixi"],
        vec!["tools", "exec", "--no-install", "--no-install", "samtools"],
        vec!["tools", "setup", "goatools", "--help"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_biov-rs"))
            .args(args)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
    }
    let output = Command::new(env!("CARGO_BIN_EXE_biov-rs"))
        .args(["tools", "--help"])
        .output()
        .unwrap();
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8(output.stderr).unwrap().contains("setup"));
}
