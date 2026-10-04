#![cfg(all(target_os = "linux", target_arch = "x86_64"))]
use biov_tools::ToolStore;
use std::{fs, path::PathBuf};
mod support;

fn fixture() -> (tempfile::TempDir, ToolStore, PathBuf, PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let pixi = support::global_manager(dir.path(), "fake-global-pixi");
    let uv = support::global_manager(dir.path(), "fake-global-uv");
    let store = ToolStore::new(dir.path().join("environments"), pixi).unwrap();
    let python = dir.path().join("paired-python");
    fs::write(&python, b"fixture interpreter").unwrap();
    (dir, store, uv, python)
}

#[test]
fn empty_list_does_not_require_managers_or_create_roots() {
    let (dir, store, uv, _) = fixture();
    fs::remove_file(dir.path().join("fake-global-pixi")).unwrap();
    fs::remove_file(&uv).unwrap();
    assert!(store
        .list_global(Some(&uv))
        .unwrap()
        .contains("No native tools installed"));
    assert!(!dir.path().join("environments").exists());
}

#[test]
fn locked_exec_setup_is_not_a_global_install() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path().join("environments");
    let locked = ToolStore::new(root.clone(), support::manager(dir.path())).unwrap();
    locked.setup("samtools").unwrap();
    fs::remove_file(dir.path().join("fake-pixi")).unwrap();
    assert!(locked
        .list_global(None)
        .unwrap()
        .contains("No native tools installed"));
    assert!(locked.workspace().is_dir());
    assert!(!root.join("pixi-global").exists());
}

#[test]
fn samtools_uses_native_manifest_and_backend_commands_without_registry_or_runner() {
    let (dir, store, uv, _) = fixture();
    store.install_global("samtools", None, None).unwrap();
    let root = dir.path().join("environments");
    assert!(root
        .join("pixi-global/manifests/pixi-global.toml")
        .is_file());
    assert!(store.bin_paths()[0].join("samtools").is_file());
    assert!(root.join("pixi-global/envs/samtools").is_dir());
    assert!(!root.join("installed").exists());
    assert!(!root.join("runners").exists());
    assert!(!store.workspace().exists());
    let listed = store.list_global(Some(&uv)).unwrap();
    assert!(listed.contains("Pixi global (") && listed.contains("samtools 1.24"));
    assert!(!listed.contains("uv tool ("));
    store.install_global("samtools", None, None).unwrap();
}

#[test]
fn goatools_delegates_all_pinned_native_entrypoints_and_isolated_roots_to_uv() {
    let (dir, store, uv, python) = fixture();
    store
        .install_global("goatools", Some(&python), Some(&uv))
        .unwrap();
    let root = dir.path().join("environments/uv-tools");
    for name in [
        "goatools",
        "find_enrichment.py",
        "go_plot.py",
        "map_to_slim.py",
        "ncbi_gene_results_to_python.py",
        "plot_go_term.py",
        "wr_hier.py",
    ] {
        let command = root.join("bin").join(name);
        assert!(fs::symlink_metadata(&command)
            .unwrap()
            .file_type()
            .is_symlink());
        assert_eq!(
            fs::canonicalize(command).unwrap(),
            root.join("tools/goatools/bin").join(name)
        );
    }
    let listed = store.list_global(Some(&uv)).unwrap();
    assert!(listed.contains("uv tool (") && listed.contains("statsmodels==0.14.6"));
    assert!(!listed.contains("Pixi global ("));
    let argv = fs::read_to_string(dir.path().join("global-argv.jsonl")).unwrap();
    let install: serde_json::Value = argv
        .lines()
        .map(|l| serde_json::from_str::<serde_json::Value>(l).unwrap())
        .find(|v| v["args"][1] == "install")
        .unwrap();
    assert_eq!(
        install["env"]["UV_TOOL_DIR"],
        root.join("tools").to_str().unwrap()
    );
    assert_eq!(
        install["env"]["UV_TOOL_BIN_DIR"],
        root.join("bin").to_str().unwrap()
    );
    assert_eq!(
        install["env"]["UV_CACHE_DIR"],
        root.join("cache").to_str().unwrap()
    );
    assert_eq!(
        install["env"]["UV_PYTHON_INSTALL_DIR"],
        root.join("python").to_str().unwrap()
    );
    assert!(argv.contains("goatools==1.6.5") && argv.contains("statsmodels==0.14.6"));
    assert!(argv.contains("--no-python-downloads") && argv.contains("--no-config"));
    store
        .install_global("goatools", Some(&python), Some(&uv))
        .unwrap();
}

#[test]
fn backend_uninstall_removes_prefix_and_commands_but_retains_scientific_workspace_cache_and_data() {
    let (dir, store, uv, python) = fixture();
    let root = dir.path().join("environments");
    let locked = ToolStore::new(root.clone(), support::manager(dir.path())).unwrap();
    locked.setup("samtools").unwrap();
    store.install_global("samtools", None, None).unwrap();
    store
        .install_global("goatools", Some(&python), Some(&uv))
        .unwrap();
    let retained = [
        root.join("pixi-global/cache/keep"),
        root.join("uv-tools/cache/keep"),
        dir.path().join("biological-data"),
        dir.path().join("analysis-output"),
    ];
    for path in &retained {
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, b"retained").unwrap();
    }
    store.uninstall_global("samtools", None).unwrap();
    store.uninstall_global("goatools", Some(&uv)).unwrap();
    assert!(!root.join("pixi-global/bin/samtools").exists());
    assert!(!root.join("pixi-global/envs/samtools").exists());
    assert!(!root.join("uv-tools/bin/goatools").exists());
    assert!(!root.join("uv-tools/tools/goatools").exists());
    assert!(locked.workspace().is_dir());
    for path in retained {
        assert_eq!(fs::read(path).unwrap(), b"retained");
    }
}

#[test]
fn foreign_samtools_command_is_refused_before_native_manifest_publication() {
    let (_dir, store, _, _) = fixture();
    let command = store.bin_paths()[0].join("samtools");
    fs::create_dir_all(command.parent().unwrap()).unwrap();
    fs::write(&command, b"foreign").unwrap();
    assert!(store.install_global("samtools", None, None).is_err());
    assert_eq!(fs::read(&command).unwrap(), b"foreign");
    assert!(!command
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("manifests/pixi-global.toml")
        .exists());
}

#[test]
fn changed_samtools_command_or_foreign_symlink_is_not_replaced_or_uninstalled() {
    let (dir, store, _, _) = fixture();
    store.install_global("samtools", None, None).unwrap();
    let command = store.bin_paths()[0].join("samtools");
    fs::write(&command, b"modified").unwrap();
    assert!(store
        .install_global("samtools", None, None)
        .unwrap_err()
        .contains("foreign"));
    assert!(store
        .uninstall_global("samtools", None)
        .unwrap_err()
        .contains("foreign"));
    assert_eq!(fs::read(&command).unwrap(), b"modified");
    fs::remove_file(&command).unwrap();
    let foreign = dir.path().join("foreign");
    fs::write(&foreign, b"foreign").unwrap();
    std::os::unix::fs::symlink(&foreign, &command).unwrap();
    assert!(store.uninstall_global("samtools", None).is_err());
    assert_eq!(fs::read(foreign).unwrap(), b"foreign");
}

#[test]
fn every_goatools_entrypoint_has_collision_preflight() {
    let (_dir, store, uv, python) = fixture();
    let target = store.bin_paths()[1].join("find_enrichment.py");
    fs::create_dir_all(target.parent().unwrap()).unwrap();
    fs::write(&target, b"foreign enrichment command").unwrap();
    assert!(store
        .install_global("goatools", Some(&python), Some(&uv))
        .unwrap_err()
        .contains("foreign"));
    assert_eq!(fs::read(&target).unwrap(), b"foreign enrichment command");
    assert!(!target
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("tools")
        .exists());
}

#[test]
fn uv_regular_or_foreign_symlink_commands_prevent_uninstall() {
    let (dir, store, uv, python) = fixture();
    store
        .install_global("goatools", Some(&python), Some(&uv))
        .unwrap();
    let target = store.bin_paths()[1].join("goatools");
    fs::remove_file(&target).unwrap();
    let foreign = dir.path().join("foreign");
    fs::write(&foreign, b"foreign").unwrap();
    std::os::unix::fs::symlink(&foreign, &target).unwrap();
    assert!(store
        .uninstall_global("goatools", Some(&uv))
        .unwrap_err()
        .contains("foreign"));
    assert_eq!(fs::read(&foreign).unwrap(), b"foreign");
    assert!(dir
        .path()
        .join("environments/uv-tools/tools/goatools")
        .is_dir());
}

#[test]
fn goatools_requires_explicit_paired_interpreter_without_ambient_python_fallback() {
    let (dir, store, uv, _) = fixture();
    assert!(store
        .install_global("goatools", None, Some(&uv))
        .unwrap_err()
        .contains("paired"));
    assert!(store
        .install_global(
            "goatools",
            Some(&dir.path().join("missing-python")),
            Some(&uv)
        )
        .is_err());
    assert!(!dir.path().join("environments").exists());
    assert!(!dir.path().join("global-argv.jsonl").exists());
}

#[test]
fn redirected_backend_directory_and_dangling_symlink_are_rejected() {
    let (dir, store, uv, python) = fixture();
    fs::create_dir(dir.path().join("environments")).unwrap();
    let foreign = dir.path().join("foreign-root");
    fs::create_dir(&foreign).unwrap();
    std::os::unix::fs::symlink(&foreign, dir.path().join("environments/pixi-global")).unwrap();
    assert!(store
        .install_global("samtools", None, None)
        .unwrap_err()
        .contains("redirected"));
    std::os::unix::fs::symlink(
        dir.path().join("missing-root"),
        dir.path().join("environments/uv-tools"),
    )
    .unwrap();
    assert!(store
        .install_global("goatools", Some(&python), Some(&uv))
        .unwrap_err()
        .contains("redirected"));
    assert_eq!(fs::read_dir(foreign).unwrap().count(), 0);
}

#[test]
fn upstream_failure_is_reported_without_an_installation_registry() {
    let (dir, store, uv, python) = fixture();
    fs::write(dir.path().join("global-fail"), b"fail").unwrap();
    assert!(store
        .install_global("samtools", None, None)
        .unwrap_err()
        .contains("failed"));
    assert!(store
        .install_global("goatools", Some(&python), Some(&uv))
        .unwrap_err()
        .contains("failed"));
    assert!(!dir.path().join("environments/installed").exists());
    assert!(!dir.path().join("environments/runners").exists());
}

#[test]
fn installed_inventory_requires_its_backend_but_not_other_manager() {
    let (dir, store, uv, _) = fixture();
    store.install_global("samtools", None, None).unwrap();
    fs::remove_file(&uv).unwrap();
    assert!(store.list_global(Some(&uv)).unwrap().contains("samtools"));
    fs::remove_file(dir.path().join("fake-global-pixi")).unwrap();
    assert!(store.list_global(Some(&uv)).is_err());
}

#[test]
fn supported_global_names_and_backend_paths_are_bounded() {
    let (dir, store, uv, python) = fixture();
    assert_eq!(
        store.bin_paths(),
        vec![
            dir.path().join("environments/pixi-global/bin"),
            dir.path().join("environments/uv-tools/bin")
        ]
    );
    assert!(store
        .install_global("../foreign", Some(&python), Some(&uv))
        .is_err());
    assert!(store.uninstall_global("anything", Some(&uv)).is_err());
    assert!(!dir.path().join("environments").exists());
}
