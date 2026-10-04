#![cfg(all(target_os = "linux", target_arch = "x86_64"))]
use biov_tools::{workspace_identity, ToolStore};
use fs4::fs_std::FileExt;
use std::{ffi::OsString, fs};
mod support;

struct Fixture {
    dir: tempfile::TempDir,
    store: ToolStore,
}
impl Fixture {
    fn new() -> Self {
        let dir = tempfile::tempdir().unwrap();
        let manager = support::manager(dir.path());
        let store = ToolStore::new(dir.path().join("environments"), manager).unwrap();
        Self { dir, store }
    }
    fn log(&self) -> String {
        fs::read_to_string(self.dir.path().join("argv.jsonl")).unwrap_or_default()
    }
}

#[test]
fn exact_workspace_identity_agrees_with_python_contract() {
    assert_eq!(
        workspace_identity(),
        "1e2ef2617576596b94e66604279be7a1fb877a1d805ae6c00d9c66fe8501f7f6"
    );
}
#[test]
fn setup_reuses_locked_install_across_instances_and_passes_literal_arguments() {
    let f = Fixture::new();
    assert_eq!(f.store.inspect("goatools").unwrap().status, "unavailable");
    let receipt = f.store.setup("goatools").unwrap();
    let same = f.store.setup("goatools").unwrap();
    assert_eq!(receipt, same);
    assert_eq!(
        f.log()
            .lines()
            .filter(|line| line.starts_with("[\"install\""))
            .count(),
        1
    );
    let arguments: Vec<OsString> = [
        "",
        "a b",
        "'quote'",
        "$(touch must-not-exist)",
        "--environment",
        "different",
        "--",
        "β",
    ]
    .into_iter()
    .map(Into::into)
    .collect();
    let status = f
        .store
        .execute("goatools", &arguments, Some(f.dir.path()))
        .unwrap();
    assert_eq!(status.code(), Some(37));
    let actual: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(f.dir.path().join("native.json")).unwrap())
            .unwrap();
    assert_eq!(
        actual["argv"],
        serde_json::json!([
            receipt.prefix.join("bin/goatools").to_str().unwrap(),
            "",
            "a b",
            "'quote'",
            "$(touch must-not-exist)",
            "--environment",
            "different",
            "--",
            "β"
        ])
    );
    assert_eq!(actual["cwd"], f.dir.path().to_str().unwrap());
    assert!(!f.dir.path().join("must-not-exist").exists());
    assert_eq!(
        f.store.inspect("goatools").unwrap().status,
        "setup_recorded"
    );
}
#[test]
fn no_install_cannot_create_workspace_or_run_unrecorded_prefix() {
    let f = Fixture::new();
    assert!(f.store.execute("samtools", &[], None).is_err());
    assert!(!f.store.workspace().exists());
    assert!(!f.log().contains("\"install\""));
    assert!(!f.log().contains("\"run\""));
}
#[test]
fn changed_manager_lock_receipt_or_prefix_never_becomes_a_ready_hit() {
    let f = Fixture::new();
    f.store.setup("samtools").unwrap();
    fs::write(f.dir.path().join("version"), "pixi 0.82.0").unwrap();
    assert!(f
        .store
        .execute("samtools", &[], None)
        .unwrap_err()
        .contains("expected Pixi"));
    fs::remove_file(f.dir.path().join("version")).unwrap();
    let marker = f
        .store
        .workspace()
        .join(".pixi/envs/samtools/conda-meta/pixi");
    let original = fs::read(&marker).unwrap();
    fs::write(&marker, b"{}").unwrap();
    assert_eq!(f.store.inspect("samtools").unwrap().status, "unavailable");
    assert!(f.store.execute("samtools", &[], None).is_err());
    fs::write(&marker, original).unwrap();
    fs::write(f.store.workspace().join("pixi.lock"), "modified").unwrap();
    assert!(f.store.execute("samtools", &[], None).is_err());
    assert!(f.store.setup("samtools").is_err());
    assert!(!f.log().contains("\"run\""));
}
#[test]
fn failed_setup_and_busy_workspace_preserve_data_without_ready_claims() {
    let f = Fixture::new();
    fs::write(f.dir.path().join("fail-install"), "yes").unwrap();
    assert!(f
        .store
        .setup("goatools")
        .unwrap_err()
        .contains("setup failed"));
    assert_eq!(f.store.inspect("goatools").unwrap().status, "unavailable");
    fs::remove_file(f.dir.path().join("fail-install")).unwrap();
    f.store.setup("goatools").unwrap();
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .open(
            f.dir
                .path()
                .join("environments/native-tool-locks")
                .join(format!("{}.lock", workspace_identity())),
        )
        .unwrap();
    FileExt::lock_exclusive(&lock).unwrap();
    assert!(f
        .store
        .execute("goatools", &[], None)
        .unwrap_err()
        .contains("busy"));
    assert!(f.store.setup("goatools").unwrap_err().contains("busy"));
    drop(lock);
    assert_eq!(
        f.store.inspect("goatools").unwrap().status,
        "setup_recorded"
    );
}
#[test]
fn unsupported_tools_fail_before_manager_or_filesystem_changes() {
    let f = Fixture::new();
    for name in ["diffdock", "conda:samtools", "../../outside", ""] {
        assert!(f.store.setup(name).is_err());
        assert!(f.store.inspect(name).is_err());
    }
    assert!(f.log().is_empty());
    assert!(!f.store.workspace().exists());
}

#[test]
fn normalized_root_and_symlinked_parent_match_manager_prefixes() {
    let f = Fixture::new();
    std::os::unix::fs::symlink(f.dir.path(), f.dir.path().join("alias")).unwrap();
    let through_alias = ToolStore::new(
        f.dir.path().join("alias/./environments"),
        f.dir.path().join("fake-pixi"),
    )
    .unwrap();
    assert_eq!(through_alias.workspace(), f.store.workspace());
    through_alias.setup("goatools").unwrap();
    assert_eq!(
        f.store.inspect("goatools").unwrap().status,
        "setup_recorded"
    );
}

#[test]
fn shared_run_lock_allows_reuse_but_missing_entrypoint_prevents_host_fallback() {
    let f = Fixture::new();
    let receipt = f.store.setup("samtools").unwrap();
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .open(
            f.dir
                .path()
                .join("environments/native-tool-locks")
                .join(format!("{}.lock", workspace_identity())),
        )
        .unwrap();
    FileExt::lock_shared(&lock).unwrap();
    assert_eq!(f.store.setup("samtools").unwrap(), receipt);
    assert_eq!(
        f.store.execute("samtools", &[], None).unwrap().code(),
        Some(37)
    );
    fs::remove_file(receipt.prefix.join("bin/samtools")).unwrap();
    assert!(f
        .store
        .execute("samtools", &[], None)
        .unwrap_err()
        .contains("executable is unavailable"));
    assert_eq!(f.store.inspect("samtools").unwrap().status, "unavailable");
}
#[test]
fn non_utf8_arguments_fail_clearly_before_manager_execution() {
    use std::os::unix::ffi::OsStringExt;
    let f = Fixture::new();
    assert!(f
        .store
        .execute("samtools", &[OsString::from_vec(vec![0xff])], None)
        .unwrap_err()
        .contains("UTF-8"));
    assert!(f.log().is_empty());
}

#[test]
fn redirected_manager_prefix_fails_before_installation() {
    let f = Fixture::new();
    fs::write(f.dir.path().join("redirect"), "yes").unwrap();
    let error = f.store.setup("goatools").unwrap_err();
    assert!(error.contains("outside"), "{error}");
    assert!(!f.log().contains("\"install\""));
    assert!(!f.dir.path().join("external").exists());
}

#[test]
fn cached_execution_rejects_later_prefix_redirect_and_symlink() {
    let f = Fixture::new();
    let receipt = f.store.setup("samtools").unwrap();
    fs::write(f.dir.path().join("redirect"), "yes").unwrap();
    assert!(f
        .store
        .execute("samtools", &[], None)
        .unwrap_err()
        .contains("outside"));
    assert!(f.store.setup("samtools").unwrap_err().contains("outside"));
    assert!(!f.log().contains("\"run\""));
    fs::remove_file(f.dir.path().join("redirect")).unwrap();
    let outside = f.dir.path().join("outside-prefix");
    fs::rename(&receipt.prefix, &outside).unwrap();
    std::os::unix::fs::symlink(&outside, &receipt.prefix).unwrap();
    assert!(f.store.execute("samtools", &[], None).is_err());
    assert_eq!(f.store.inspect("samtools").unwrap().status, "unavailable");
    assert!(!f.log().contains("\"run\""));
}
