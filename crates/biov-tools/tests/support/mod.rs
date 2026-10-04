//! Executable protocol fixtures independent of source-distribution file modes.
use std::{
    fs,
    os::unix::fs::PermissionsExt,
    path::{Path, PathBuf},
    sync::OnceLock,
};

static FIXTURE: OnceLock<tempfile::TempDir> = OnceLock::new();

/// Complete the one script write before any fixture subprocess can start.
/// Parallel child creation can inherit another thread's live write descriptor,
/// producing transient ETXTBSY even after the originating writer closes it.
pub fn prepare() -> &'static tempfile::TempDir {
    FIXTURE.get_or_init(|| {
        let directory = tempfile::tempdir().unwrap();
        let executable = directory.path().join("fixture.py");
        fs::write(&executable, include_bytes!("../fixtures/fake_pixi.py")).unwrap();
        fs::set_permissions(executable, fs::Permissions::from_mode(0o755)).unwrap();
        directory
    })
}

/// Each symlink retains its own argv[0]/__file__ directory for isolated records.
pub fn manager(directory: &Path) -> PathBuf {
    let executable = prepare().path().join("fixture.py");
    let manager = directory.join("fake-pixi");
    std::os::unix::fs::symlink(executable, &manager).unwrap();
    manager
}

#[allow(dead_code)] // This module is also included by locked-only integration suites.
static GLOBAL_FIXTURE: OnceLock<tempfile::TempDir> = OnceLock::new();
/// Upstream global/tool CLI fixture, kept separate from locked scientific exec.
#[allow(dead_code)] // Used by the global backend suite, not locked-only suites.
pub fn global_manager(directory: &Path, name: &str) -> PathBuf {
    let fixture = GLOBAL_FIXTURE.get_or_init(|| {
        let directory = tempfile::tempdir().unwrap();
        let executable = directory.path().join("global.py");
        fs::write(&executable, include_bytes!("../fixtures/fake_global.py")).unwrap();
        fs::set_permissions(executable, fs::Permissions::from_mode(0o755)).unwrap();
        directory
    });
    let manager = directory.join(name);
    std::os::unix::fs::symlink(fixture.path().join("global.py"), &manager).unwrap();
    manager
}
