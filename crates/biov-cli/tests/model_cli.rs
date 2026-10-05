#![cfg(all(target_os = "linux", target_arch = "x86_64"))]
use std::{
    fs,
    os::unix::fs::PermissionsExt,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::OnceLock,
};

const REVISION: &str = "f171d7baecaf37b5da5a3616d8833b9969753535";
static CLIENT: OnceLock<tempfile::TempDir> = OnceLock::new();
fn client() -> PathBuf {
    CLIENT
        .get_or_init(|| {
            let directory = tempfile::tempdir().unwrap();
            let file = directory.path().join("fixture.py");
            fs::write(&file, include_bytes!("fixtures/fake_hf.py")).unwrap();
            fs::set_permissions(&file, fs::Permissions::from_mode(0o755)).unwrap();
            directory
        })
        .path()
        .join("fixture.py")
}
struct Fixture {
    directory: tempfile::TempDir,
    hf: PathBuf,
    output: PathBuf,
}
impl Fixture {
    fn new() -> Self {
        let directory = tempfile::tempdir().unwrap();
        let hf = directory.path().join("hf");
        std::os::unix::fs::symlink(client(), &hf).unwrap();
        let output = directory.path().join("model O'Connor with spaces");
        Self {
            directory,
            hf,
            output,
        }
    }
    fn command(&self) -> Command {
        let mut command = Command::new(env!("CARGO_BIN_EXE_biov"));
        command
            .current_dir(self.directory.path())
            .env("BIOV_FAKE_CLIENT_HOME", self.directory.path())
            .env_remove("BIOV_HF_BIN")
            .env_remove("BIOV_UV_BIN")
            .env_remove("HF_ENDPOINT")
            .env_remove("HUGGINGFACE_CO_STAGING");
        command
    }
    fn download(&self, files: &[&str]) -> Output {
        self.command()
            .args(["model", "download", "--revision", REVISION, "--local-dir"])
            .arg(&self.output)
            .arg("--hf")
            .arg(&self.hf)
            .arg("example/model")
            .args(files)
            .output()
            .unwrap()
    }
    fn calls(&self) -> String {
        fs::read_to_string(self.directory.path().join("calls.jsonl")).unwrap_or_default()
    }
}
fn result(output: &Output) -> serde_json::Value {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}

#[test]
fn exact_download_portable_offline_reuse_and_damage_rejection() {
    let f = Fixture::new();
    let files = [
        "config.json",
        "nested/tokenizer O'Connor.json",
        "$(touch never-executed).json",
    ];
    let output = f.download(&files);
    let record = result(&output);
    assert_eq!(record["status"], "downloaded");
    assert_eq!(record["resource"]["repository"], "example/model");
    assert_eq!(record["resource"]["revision"], REVISION);
    assert_eq!(record["resource"]["inventory"].as_array().unwrap().len(), 3);
    let diagnostics = String::from_utf8(output.stderr).unwrap();
    assert!(
        diagnostics.contains("upstream download output")
            && diagnostics.contains("upstream progress")
    );
    assert!(!f.directory.path().join("never-executed").exists());
    let calls = f.calls();
    assert_eq!(calls.lines().count(), 3);
    let repeat = result(&f.download(&files));
    assert_eq!(repeat["status"], "reused");
    assert_eq!(f.calls(), calls);
    let moved = f.directory.path().join("moved resource");
    fs::rename(&f.output, &moved).unwrap();
    let offline = f
        .command()
        .args(["model", "inspect"])
        .arg(&moved)
        .env("PATH", "")
        .env("BIOV_HF_BIN", "/missing/hf")
        .env("BIOV_UV_BIN", "/missing/uv")
        .output()
        .unwrap();
    assert_eq!(result(&offline)["status"], "verified");
    assert_eq!(f.calls(), calls);
    fs::write(moved.join("config.json"), b"modified").unwrap();
    let damaged = f
        .command()
        .args(["model", "inspect"])
        .arg(&moved)
        .output()
        .unwrap();
    assert_eq!(damaged.status.code(), Some(2));
    assert!(damaged.stdout.is_empty());
}

#[test]
fn missing_selected_file_and_backend_failure_never_publish() {
    for backend_failure in [false, true] {
        let f = Fixture::new();
        fs::write(
            f.directory.path().join(if backend_failure {
                "fail-download"
            } else {
                "omit-file"
            }),
            "missing.json",
        )
        .unwrap();
        let output = f.download(&["config.json", "missing.json"]);
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
        assert!(!f.output.exists());
        let error = String::from_utf8(output.stderr).unwrap();
        assert!(error.contains("unpublished files retained"));
        if backend_failure {
            assert!(error.contains("23"));
        }
        let retained: Vec<_> = fs::read_dir(f.directory.path())
            .unwrap()
            .map(|e| e.unwrap().path())
            .filter(|p| {
                p.file_name()
                    .unwrap()
                    .to_string_lossy()
                    .starts_with(".biov-model-download-")
            })
            .collect();
        assert_eq!(retained.len(), 1);
        assert!(retained[0]
            .join("BIOV_INCOMPLETE_MODEL_DOWNLOAD.json")
            .is_file());
        assert!(!retained[0].join("BIOV_MODEL_RESOURCE.json").exists());
    }
}

#[test]
fn invalid_selection_and_existing_destination_have_no_backend_side_effects() {
    let f = Fixture::new();
    for files in [
        vec![],
        vec!["../escape"],
        vec!["*.bin"],
        vec!["folder/"],
        vec!["config.json", "config.json"],
    ] {
        let output = f.download(&files);
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
    }
    assert!(f.calls().is_empty());
    assert!(!f.output.exists());
    fs::create_dir(&f.output).unwrap();
    fs::write(f.output.join("keep.json"), b"external data").unwrap();
    let output = f.download(&["config.json"]);
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(
        fs::read(f.output.join("keep.json")).unwrap(),
        b"external data"
    );
    assert!(f.calls().is_empty());
}

#[test]
fn official_uv_fallback_and_relative_installed_executable() {
    let f = Fixture::new();
    let uv = f.directory.path().join("uv");
    std::os::unix::fs::symlink(client(), &uv).unwrap();
    let python = Command::new("python3")
        .args(["-c", "import sys; print(sys.executable)"])
        .output()
        .unwrap();
    let python = Path::new(std::str::from_utf8(&python.stdout).unwrap().trim());
    let bin = f.directory.path().join("only-uv");
    fs::create_dir(&bin).unwrap();
    std::os::unix::fs::symlink(&uv, bin.join("uv")).unwrap();
    std::os::unix::fs::symlink(python, bin.join("python3")).unwrap();
    let output = f
        .command()
        .args(["model", "download", "--revision", REVISION, "--local-dir"])
        .arg(&f.output)
        .arg("example/model")
        .arg("config.json")
        .env("PATH", &bin)
        .output()
        .unwrap();
    assert_eq!(result(&output)["status"], "downloaded");
    let lines: Vec<serde_json::Value> = f
        .calls()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(lines.len(), 3);
    assert!(lines.iter().all(|entry| entry["backend"] == "uv"));
    let another = Fixture::new();
    let output = another
        .command()
        .current_dir(another.directory.path())
        .args(["model", "download", "--revision", REVISION, "--local-dir"])
        .arg(&another.output)
        .args(["--hf", "./hf", "example/model", "config.json"])
        .output()
        .unwrap();
    assert_eq!(result(&output)["status"], "downloaded");
}

#[test]
fn no_install_revision_conflicts_and_help_are_explicit() {
    let f = Fixture::new();
    let _ = result(&f.download(&["config.json"]));
    let before = f.calls();
    let conflict = f.download(&["other.json"]);
    assert_eq!(conflict.status.code(), Some(2));
    assert_eq!(f.calls(), before);
    for words in [
        vec!["model", "--help"],
        vec!["model", "download", "--help"],
        vec!["model", "inspect", "--help"],
    ] {
        let output = f.command().args(words).env("PATH", "").output().unwrap();
        assert!(output.status.success());
        assert!(output.stdout.is_empty());
    }
    let output = f
        .command()
        .args([
            "model",
            "download",
            "--no-install",
            "--revision",
            REVISION,
            "--local-dir",
        ])
        .arg(f.directory.path().join("missing-client"))
        .args(["example/model", "config.json"])
        .env("PATH", "")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    assert_eq!(f.calls(), before);
}

#[test]
fn relative_destinations_and_source_endpoint_are_preflighted() {
    for directory in ["model", "models/nested resource"] {
        let f = Fixture::new();
        let output = f
            .command()
            .current_dir(f.directory.path())
            .args([
                "model",
                "download",
                "--revision",
                REVISION,
                "--local-dir",
                directory,
                "--hf",
                "./hf",
                "example/model",
                "config.json",
            ])
            .output()
            .unwrap();
        assert_eq!(result(&output)["status"], "downloaded");
        assert!(f
            .directory
            .path()
            .join(directory)
            .join("config.json")
            .is_file());
    }
    let f = Fixture::new();
    let output = f
        .command()
        .args(["model", "download", "--revision", REVISION, "--local-dir"])
        .arg(&f.output)
        .arg("--hf")
        .arg(&f.hf)
        .args(["example/model", "config.json"])
        .env("HF_ENDPOINT", "https://different.example.invalid")
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("custom HF_ENDPOINT"));
    assert!(f.calls().is_empty());
    assert!(!f.output.exists());
}

#[test]
fn concurrent_destination_is_preserved_without_replacing_or_adopting_it() {
    let f = Fixture::new();
    fs::write(
        f.directory.path().join("concurrent-destination"),
        f.output.to_str().unwrap(),
    )
    .unwrap();
    let output = f.download(&["config.json"]);
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    assert_eq!(
        fs::read(f.output.join("keep.txt")).unwrap(),
        b"concurrent external data"
    );
    assert!(!f.output.join("config.json").exists());
    assert!(!f.output.join("BIOV_MODEL_RESOURCE.json").exists());
    let retained: Vec<_> = fs::read_dir(f.directory.path())
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| {
            p.file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with(".biov-model-download-")
        })
        .collect();
    assert_eq!(retained.len(), 1);
    let inspect = f
        .command()
        .args(["model", "inspect"])
        .arg(&retained[0])
        .env("PATH", "")
        .output()
        .unwrap();
    assert_eq!(result(&inspect)["status"], "verified");
}

#[test]
fn symlinked_destination_ancestors_are_rejected_before_backend_access() {
    let f = Fixture::new();
    let real = f.directory.path().join("real parent");
    let alias = f.directory.path().join("alias parent");
    fs::create_dir(&real).unwrap();
    std::os::unix::fs::symlink(&real, &alias).unwrap();
    for destination in [alias.join("model"), alias.join("nested/model")] {
        let output = f
            .command()
            .args(["model", "download", "--revision", REVISION, "--local-dir"])
            .arg(&destination)
            .arg("--hf")
            .arg(&f.hf)
            .args(["example/model", "config.json"])
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8(output.stderr)
            .unwrap()
            .contains("real directories"));
        assert!(!destination.exists());
        assert!(f.calls().is_empty());
    }
}

#[test]
fn staging_source_override_is_rejected_and_disabled_mode_is_allowed() {
    for value in ["1", "ON", "yes", "true"] {
        let f = Fixture::new();
        let output = f
            .command()
            .args(["model", "download", "--revision", REVISION, "--local-dir"])
            .arg(&f.output)
            .arg("--hf")
            .arg(&f.hf)
            .args(["example/model", "config.json"])
            .env("HUGGINGFACE_CO_STAGING", value)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8(output.stderr)
            .unwrap()
            .contains("HUGGINGFACE_CO_STAGING"));
        assert!(!f.output.exists());
        assert!(f.calls().is_empty());
    }
    let f = Fixture::new();
    let output = f
        .command()
        .args(["model", "download", "--revision", REVISION, "--local-dir"])
        .arg(&f.output)
        .arg("--hf")
        .arg(&f.hf)
        .args(["example/model", "config.json"])
        .env("HUGGINGFACE_CO_STAGING", "0")
        .output()
        .unwrap();
    assert_eq!(result(&output)["status"], "downloaded");
}

#[test]
fn installed_client_precedence_and_incompatible_client_policy() {
    for (version, explicit, no_install, expected_backends) in [
        ("0.34.0", false, false, vec!["hf", "hf", "hf"]),
        ("0.33.9", false, false, vec!["hf", "uv", "uv", "uv"]),
        ("0.33.9", true, false, vec!["hf"]),
        ("0.33.9", false, true, vec!["hf"]),
    ] {
        let f = Fixture::new();
        fs::write(f.directory.path().join("version"), version).unwrap();
        let bin = f.directory.path().join("client bin");
        fs::create_dir(&bin).unwrap();
        for name in ["hf", "uv"] {
            std::os::unix::fs::symlink(client(), bin.join(name)).unwrap();
        }
        let python = Command::new("python3")
            .args(["-c", "import sys; print(sys.executable)"])
            .output()
            .unwrap();
        std::os::unix::fs::symlink(
            Path::new(std::str::from_utf8(&python.stdout).unwrap().trim()),
            bin.join("python3"),
        )
        .unwrap();
        let mut command = f.command();
        command.args(["model", "download", "--revision", REVISION, "--local-dir"]);
        command.arg(&f.output);
        if explicit {
            command.arg("--hf").arg(&f.hf);
        }
        if no_install {
            command.arg("--no-install");
        }
        let output = command
            .args(["example/model", "config.json"])
            .env("PATH", &bin)
            .output()
            .unwrap();
        if explicit || no_install {
            assert_eq!(output.status.code(), Some(2));
            assert!(output.stdout.is_empty());
            assert!(!f.output.exists());
        } else {
            assert_eq!(result(&output)["status"], "downloaded");
        }
        let backends: Vec<String> = f
            .calls()
            .lines()
            .map(|line| {
                serde_json::from_str::<serde_json::Value>(line).unwrap()["backend"]
                    .as_str()
                    .unwrap()
                    .to_owned()
            })
            .collect();
        assert_eq!(backends, expected_backends);
    }
}

#[test]
fn failed_download_retry_uses_distinct_fresh_staging() {
    let f = Fixture::new();
    fs::write(f.directory.path().join("fail-download"), "").unwrap();
    for _ in 0..2 {
        assert_eq!(f.download(&["config.json"]).status.code(), Some(2));
    }
    let directories: Vec<String> = f
        .calls()
        .lines()
        .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
        .filter(|call| call["args"][0] == "download" && call["args"][1] != "--help")
        .map(|call| call["local_dir"].as_str().unwrap().to_owned())
        .collect();
    assert_eq!(directories.len(), 2);
    assert_ne!(directories[0], directories[1]);
    assert!(directories
        .iter()
        .all(|directory| Path::new(directory).is_dir()));
    assert!(!f.output.exists());
}

#[test]
fn caller_cwd_preserves_relative_upstream_environment_paths() {
    let f = Fixture::new();
    fs::write(f.directory.path().join("relative-environment"), "").unwrap();
    let mut command = f.command();
    for name in [
        "HF_HOME",
        "HF_TOKEN_PATH",
        "HF_HUB_CACHE",
        "UV_TOOL_DIR",
        "UV_CACHE_DIR",
    ] {
        fs::write(f.directory.path().join(name), "existing caller state").unwrap();
        command.env(name, name);
    }
    let output = command
        .args(["model", "download", "--revision", REVISION, "--local-dir"])
        .arg(&f.output)
        .arg("--hf")
        .arg(&f.hf)
        .args(["example/model", "config.json"])
        .output()
        .unwrap();
    assert_eq!(result(&output)["status"], "downloaded");
    for call in f.calls().lines() {
        let call: serde_json::Value = serde_json::from_str(call).unwrap();
        assert_eq!(call["cwd"], f.directory.path().to_str().unwrap());
    }
}
