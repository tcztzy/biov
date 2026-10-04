//! Opt-in regressions against actual pinned Pixi, independent of fake-manager tests.
//! Synthetic prefix tests establish activation/argv behavior, not package installation.
use super::*;
use std::os::unix::{fs::PermissionsExt, process::ExitStatusExt};

fn manager() -> PathBuf {
    PathBuf::from(
        std::env::var_os("BIOV_TEST_REAL_PIXI").expect("set BIOV_TEST_REAL_PIXI to Pixi 0.81.0"),
    )
}
fn shell_word(path: &Path) -> String {
    format!("'{}'", path.to_str().unwrap().replace('\'', "'\"'\"'"))
}
fn synthetic_store(root: &Path) -> (ToolStore, PathBuf) {
    let store = ToolStore::new(root.to_owned(), manager()).unwrap();
    store.check_manager().unwrap();
    let _guard = store.operation_lock(true, true).unwrap();
    store.publish_workspace().unwrap();
    let prefix = store.workspace().join(".pixi/envs/samtools");
    fs::create_dir_all(prefix.join("bin")).unwrap();
    fs::create_dir_all(prefix.join("conda-meta")).unwrap();
    fs::write(
        prefix.join("conda-meta/pixi"),
        serde_json::to_vec(&serde_json::json!({
            "manifest_path": store.workspace().join("pyproject.toml"),
            "environment_name": "samtools", "pixi_version": PIXI_VERSION,
            "environment_lock_file_hash": "0000000000000000",
            "resolved_platform": {"subdir": "linux-64", "virtual_packages": []},
        "minimum_supported_platform": {"subdir": "linux-64", "requirements": []}
        }))
        .unwrap(),
    )
    .unwrap();
    let record = prefix.join("observed");
    let script = format!("#!/bin/sh\nprintf '%s\\n' \"$#\" > {}\nif [ \"$#\" -gt 0 ]; then printf '%s\\0' \"$@\"; fi > {}\nif [ \"${{1-}}\" = '__signal__' ]; then kill -TERM $$; fi\nexit 37\n", shell_word(&record.with_extension("count")), shell_word(&record.with_extension("args")));
    fs::write(prefix.join("bin/samtools"), script).unwrap();
    fs::set_permissions(
        prefix.join("bin/samtools"),
        fs::Permissions::from_mode(0o755),
    )
    .unwrap();
    let receipt = store.expected_receipt("samtools", prefix).unwrap();
    fs::write(
        store.receipt_path("samtools"),
        serde_json::to_vec(&receipt).unwrap(),
    )
    .unwrap();
    (store, record)
}

#[test]
#[ignore = "requires BIOV_TEST_REAL_PIXI=path/to/Pixi0.81.0; synthetic prefix, no downloads"]
fn real_pixi_activation_preserves_zero_and_literal_arguments() {
    let directory = tempfile::tempdir().unwrap();
    for name in ["plain", "root with spaces", "O'Connor", "space O'Connor"] {
        let (store, record) = synthetic_store(&directory.path().join(name));
        let status = store
            .execute("samtools", &[], Some(directory.path()))
            .unwrap();
        assert_eq!(status.code(), Some(37));
        assert_eq!(
            fs::read_to_string(record.with_extension("count")).unwrap(),
            "0\n"
        );
        assert!(fs::read(record.with_extension("args")).unwrap().is_empty());
        let arguments: Vec<OsString> = [
            "",
            "a b",
            "O'Connor",
            "$(touch must-not-exist)",
            "$HOME",
            "`echo text`",
            "a\\b",
            "line\nbreak",
            "--",
            "β",
        ]
        .into_iter()
        .map(Into::into)
        .collect();
        assert_eq!(
            store
                .execute("samtools", &arguments, Some(directory.path()))
                .unwrap()
                .code(),
            Some(37)
        );
        let expected: Vec<u8> = arguments
            .iter()
            .flat_map(|value| {
                value
                    .to_str()
                    .unwrap()
                    .as_bytes()
                    .iter()
                    .copied()
                    .chain([0])
            })
            .collect();
        assert_eq!(fs::read(record.with_extension("args")).unwrap(), expected);
        assert_eq!(
            store
                .execute("samtools", &["__signal__".into()], None)
                .unwrap()
                .signal(),
            Some(15)
        );
        assert!(!directory.path().join("must-not-exist").exists());
    }
    // Pinned Pixi can expand dollar-containing prefixes while emitting JSON.
    // BioV must fail closed instead of running under a different activation.
    let (store, record) = synthetic_store(&directory.path().join("literal-$HOME"));
    assert!(store.execute("samtools", &[], None).is_err());
    assert!(!record.with_extension("count").exists());
}

#[test]
#[ignore = "requires BIOV_TEST_REAL_PIXI=path/to/Pixi0.81.0; synthetic prefix, no downloads"]
fn real_pixi_conda_activation_exports_reach_native_launch_once() {
    let directory = tempfile::tempdir().unwrap();
    let (store, record) = synthetic_store(&directory.path().join("activation O'Connor"));
    let prefix = record.parent().unwrap();
    let activation_directory = prefix.join("etc/conda/activate.d");
    fs::create_dir_all(&activation_directory).unwrap();
    let counter = record.with_extension("activation-count");
    fs::write(
        activation_directory.join("biov-export-probe.sh"),
        format!(
            "export BIOV_REAL_PIXI_ACTIVATION_EXPORT='activated value with spaces'\nprintf x >> {}\n",
            shell_word(&counter),
        ),
    )
    .unwrap();
    let export_record = record.with_extension("activation-export");
    let prefix_record = record.with_extension("activation-prefix");
    fs::write(
        prefix.join("bin/samtools"),
        format!(
            "#!/bin/sh\nprintf '%s' \"${{BIOV_REAL_PIXI_ACTIVATION_EXPORT-}}\" > {}\nprintf '%s' \"${{CONDA_PREFIX-}}\" > {}\nexit 37\n",
            shell_word(&export_record),
            shell_word(&prefix_record),
        ),
    )
    .unwrap();

    // JSON activation already runs package hooks. The exact native executable
    // must receive their exports without BioV sourcing the hooks a second time.
    let status = store
        .execute("samtools", &[], Some(directory.path()))
        .unwrap();
    assert_eq!(status.code(), Some(37));
    assert_eq!(
        fs::read_to_string(export_record).unwrap(),
        "activated value with spaces"
    );
    assert_eq!(
        fs::read_to_string(prefix_record).unwrap(),
        prefix.to_str().unwrap()
    );
    assert_eq!(fs::read_to_string(counter).unwrap(), "x");
}

#[test]
#[ignore = "requires BIOV_TEST_REAL_PIXI and cached/network access for locked Samtools"]
fn real_pixi_repairs_missing_locked_samtools() {
    let directory = tempfile::tempdir().unwrap();
    let store = ToolStore::new(directory.path().join("repair O'Connor"), manager()).unwrap();
    let receipt = store.setup("samtools").unwrap();
    let lock = fs::read(store.workspace().join("pixi.lock")).unwrap();
    fs::remove_file(receipt.prefix.join("bin/samtools")).unwrap();
    assert!(store.execute("samtools", &[], None).is_err());
    assert_eq!(store.inspect("samtools").unwrap().status, "unavailable");
    let repaired = store.setup("samtools").unwrap();
    assert_eq!(repaired.prefix, receipt.prefix);
    assert_eq!(fs::read(store.workspace().join("pixi.lock")).unwrap(), lock);
    assert!(store
        .execute("samtools", &["--version".into()], None)
        .unwrap()
        .success());
}
