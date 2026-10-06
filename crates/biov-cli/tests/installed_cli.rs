#![cfg(all(target_os = "linux", target_arch = "x86_64"))]
//! Manager protocol fixtures only. Real installs and BioV-removal independence
//! are established by tests/test_unified_install.py against built wheels.
use std::{
    fs,
    os::unix::fs::PermissionsExt,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::{Arc, Mutex, Weak},
};

static FIXTURE: Mutex<Weak<tempfile::TempDir>> = Mutex::new(Weak::new());

fn prepare() -> Arc<tempfile::TempDir> {
    let mut shared = FIXTURE.lock().unwrap();
    if let Some(fixtures) = shared.upgrade() {
        return fixtures;
    }
    let fixtures = Arc::new({
        let directory = tempfile::tempdir().unwrap();
        let script = directory.path().join("manager.py");
        fs::write(&script, include_bytes!("fixtures/fake_global_managers.py")).unwrap();
        fs::set_permissions(&script, fs::Permissions::from_mode(0o755)).unwrap();
        fs::copy(env!("CARGO_BIN_EXE_biov"), directory.path().join("biov")).unwrap();
        fs::write(
            directory.path().join("python"),
            "protocol interpreter placeholder",
        )
        .unwrap();
        directory
    });
    *shared = Arc::downgrade(&fixtures);
    fixtures
}

fn manager(source: &tempfile::TempDir, directory: &Path, name: &str) -> PathBuf {
    let manager = directory.join(name);
    std::os::unix::fs::symlink(source.path().join("manager.py"), &manager).unwrap();
    manager
}

fn successful(output: &Output) {
    assert!(
        output.status.success(),
        "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn public_lifecycle_delegates_to_isolated_upstream_managers() {
    let fixtures = prepare();
    let dir = tempfile::tempdir().unwrap();
    let pixi = manager(&fixtures, dir.path(), "fake-pixi");
    let uv = manager(&fixtures, dir.path(), "fake-uv");
    let root = dir.path().join("environments 'literal'");
    let bin = dir.path().join("management-bin");
    fs::create_dir(&bin).unwrap();
    let cli = bin.join("biov");
    // Immutable template links create no live executable write descriptor after
    // the common preparation barrier. Both temporary roots share the filesystem.
    fs::hard_link(fixtures.path().join("biov"), &cli).unwrap();
    // Paired interpreter path selection is checked without launching Python;
    // actual installed-interpreter behavior belongs to the real wheel gate.
    fs::hard_link(fixtures.path().join("python"), bin.join("python")).unwrap();
    let log = dir.path().join("calls.jsonl");
    let foreign_pixi = dir.path().join("foreign-pixi");
    let foreign_uv = dir.path().join("foreign-uv");
    let run = |action: &str, name: Option<&str>| {
        let mut command = Command::new(&cli);
        command
            .arg(action)
            .arg("--environment-root")
            .arg(&root)
            .arg("--pixi")
            .arg(&pixi)
            .arg("--uv")
            .arg(&uv)
            .env("BIOV_TEST_GLOBAL_LOG", &log)
            .env("PIXI_HOME", &foreign_pixi)
            .env("UV_TOOL_DIR", &foreign_uv)
            .env("UV_TOOL_BIN_DIR", dir.path().join("foreign-bin"))
            .current_dir(dir.path());
        if let Some(name) = name {
            command.arg(name);
        }
        command.output().unwrap()
    };
    let empty = run("list", None);
    successful(&empty);
    assert!(String::from_utf8_lossy(&empty.stdout).contains("No native tools installed"));
    assert!(!root.exists());
    assert!(!log.exists());
    for name in ["samtools", "goatools"] {
        let output = run("install", Some(name));
        successful(&output);
        assert!(String::from_utf8_lossy(&output.stderr).contains("export PATH="));
    }
    assert!(root.join("pixi-global/bin/samtools").is_file());
    assert!(root.join("uv-tools/bin/goatools").is_file());
    assert!(!root.join("runners").exists());
    assert!(!root.join("installed").exists());
    assert!(!foreign_pixi.exists());
    assert!(!foreign_uv.exists());
    let listing = run("list", None);
    successful(&listing);
    let text = String::from_utf8_lossy(&listing.stdout);
    assert!(text.contains("Pixi global (") && text.contains("samtools 1.24"));
    assert!(text.contains("uv tool (") && text.contains("goatools v1.6.5"));
    for name in ["samtools", "goatools"] {
        successful(&run("uninstall", Some(name)));
    }
    assert!(!root.join("pixi-global/bin/samtools").exists());
    assert!(!root.join("pixi-global/envs/samtools").exists());
    assert!(!root.join("uv-tools/bin/goatools").exists());
    assert!(!root.join("uv-tools/tools/goatools").exists());

    let calls: Vec<serde_json::Value> = fs::read_to_string(log)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    let native = calls
        .iter()
        .find(|call| call["argv"][0] == "global" && call["argv"][1] == "install")
        .unwrap();
    assert_eq!(
        native["PIXI_HOME"],
        root.join("pixi-global").to_str().unwrap()
    );
    let arguments = native["argv"].as_array().unwrap();
    assert!(arguments.iter().any(|arg| arg == "samtools==1.24"));
    assert!(arguments.iter().any(|arg| arg == "--expose"));
    assert!(arguments.iter().any(|arg| arg == "--no-shortcuts"));
    let python = calls
        .iter()
        .find(|call| call["argv"][0] == "tool" && call["argv"][1] == "install")
        .unwrap();
    assert_eq!(
        python["UV_TOOL_DIR"],
        root.join("uv-tools/tools").to_str().unwrap()
    );
    assert_eq!(
        python["UV_TOOL_BIN_DIR"],
        root.join("uv-tools/bin").to_str().unwrap()
    );
    let arguments = python["argv"].as_array().unwrap();
    assert!(arguments.iter().any(|arg| arg == "goatools==1.6.5"));
    assert!(arguments.iter().any(|arg| arg == "statsmodels==0.14.6"));
    assert!(arguments
        .iter()
        .any(|arg| arg == bin.join("python").to_str().unwrap()));
    assert!(arguments.iter().any(|arg| arg == "--no-python-downloads"));
    assert!(arguments.iter().any(|arg| arg == "--no-config"));
}

#[test]
fn removed_custom_installer_options_and_private_runner_are_rejected() {
    let _fixtures = prepare();
    for args in [
        vec!["install", "--bin-dir", "/tmp/ignored", "samtools"],
        vec!["list", "--json"],
        vec!["installed-run", "samtools"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_biov"))
            .args(args)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
    }
}
