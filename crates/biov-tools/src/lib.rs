//! Local locked-tool bridge. Pixi owns resolution, installation and activation.
//! This first migration slice supports two preparation-free bundled tools.
use fs4::fs_std::FileExt;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    ffi::OsString,
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    process::{Command, ExitStatus, Stdio},
};

pub const PIXI_VERSION: &str = "0.81.0";
pub const SUPPORTED_TOOLS: &[&str] = &["samtools", "goatools"];
const MANIFEST: &[u8] = include_bytes!("../../../src/biov/assets/environments/pyproject.toml");
const LOCK: &[u8] = include_bytes!("../../../src/biov/assets/environments/pixi.lock");
const MAX_RECORD: u64 = 64 * 1024;

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
pub fn workspace_identity() -> String {
    let mut hash = Sha256::new();
    hash.update(MANIFEST);
    hash.update([0]);
    hash.update(LOCK);
    format!("{:x}", hash.finalize())
}

/// Configuration is explicit or taken from documented environment variables.
/// Native commands do not read Python configuration or route to SSH.
#[derive(Debug, Clone)]
pub struct ToolStore {
    root: PathBuf,
    pixi: PathBuf,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Receipt {
    pub format_version: u32,
    pub environment: String,
    pub workspace_sha256: String,
    pub lock_sha256: String,
    pub manager_version: String,
    pub platform: String,
    pub prefix: PathBuf,
    pub prefix_marker_sha256: String,
}
#[derive(Debug, Serialize)]
pub struct Inspection {
    pub environment: String,
    pub source: &'static str,
    pub ownership: &'static str,
    pub workspace: PathBuf,
    pub workspace_sha256: String,
    pub lock_sha256: String,
    pub manager: PathBuf,
    pub required_manager_version: &'static str,
    pub platform: &'static str,
    pub status: &'static str,
    pub detail: String,
    pub receipt: Option<Receipt>,
}

impl ToolStore {
    pub fn new(root: PathBuf, pixi: PathBuf) -> Result<Self, String> {
        if !cfg!(all(target_os = "linux", target_arch = "x86_64")) {
            return Err(
                "native locked tools currently support Linux x86_64 (linux-64) only".into(),
            );
        }
        if root.as_os_str().is_empty() || pixi.as_os_str().is_empty() {
            return Err("tool root and Pixi executable must be nonempty".into());
        }
        let root = if root.is_absolute() {
            root
        } else {
            std::env::current_dir()
                .map_err(|e| e.to_string())?
                .join(root)
        };
        if root.to_str().is_none() {
            return Err("tool root must be UTF-8".into());
        }
        let root = canonical_future_path(&root)?;
        // Resolve explicit relative executable paths before changing the child cwd.
        let pixi = if pixi.components().count() > 1 && !pixi.is_absolute() {
            std::env::current_dir()
                .map_err(|e| e.to_string())?
                .join(pixi)
        } else {
            pixi
        };
        Ok(Self { root, pixi })
    }
    pub fn from_environment(root: Option<PathBuf>, pixi: Option<PathBuf>) -> Result<Self, String> {
        if root.is_none()
            && std::env::var_os("BIOV_ENVIRONMENT_ROOT").is_none()
            && std::env::var_os("XDG_DATA_HOME").is_none()
            && std::env::var_os("HOME").is_none()
        {
            return Err("set --environment-root or BIOV_ENVIRONMENT_ROOT".into());
        }
        let root = root
            .or_else(|| std::env::var_os("BIOV_ENVIRONMENT_ROOT").map(PathBuf::from))
            .unwrap_or_else(|| {
                std::env::var_os("XDG_DATA_HOME")
                    .map(PathBuf::from)
                    .unwrap_or_else(|| {
                        PathBuf::from(std::env::var_os("HOME").unwrap_or_default())
                            .join(".local/share")
                    })
                    .join("biov/environments")
            });
        let pixi = pixi
            .or_else(|| std::env::var_os("BIOV_PIXI_BIN").map(PathBuf::from))
            .unwrap_or_else(|| {
                let managed = root.join(format!("pixi-{PIXI_VERSION}/pixi"));
                if managed.is_file() {
                    managed
                } else {
                    PathBuf::from("pixi")
                }
            });
        Self::new(root, pixi)
    }
    pub fn workspace(&self) -> PathBuf {
        self.root.join("workspaces").join(workspace_identity())
    }
    fn supported(name: &str) -> Result<(), String> {
        if SUPPORTED_TOOLS.contains(&name) {
            Ok(())
        } else {
            Err(format!("unsupported native tool {name:?}; initial tools are samtools and goatools; use the Python biov route for other environments"))
        }
    }
    fn check_manager(&self) -> Result<(), String> {
        let output = self.capture([OsString::from("--version")])?;
        if output.trim() != format!("pixi {PIXI_VERSION}") {
            return Err(format!("expected Pixi {PIXI_VERSION}, received {:?}; set BIOV_PIXI_BIN to the matching existing executable", output.trim()));
        }
        Ok(())
    }
    fn capture(&self, args: impl IntoIterator<Item = OsString>) -> Result<String, String> {
        let mut child = Command::new(&self.pixi)
            .args(args)
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .map_err(|e| {
                format!(
                    "cannot run Pixi {:?}: {e}; install Pixi {PIXI_VERSION} or set BIOV_PIXI_BIN",
                    self.pixi
                )
            })?;
        let mut bytes = Vec::new();
        let result = child
            .stdout
            .take()
            .ok_or("missing Pixi stdout")?
            .take(1024 * 1024 + 1)
            .read_to_end(&mut bytes);
        if result.is_err() || bytes.len() > 1024 * 1024 {
            let _ = child.kill();
            let _ = child.wait();
            return Err("Pixi inspection output exceeds 1 MiB or could not be read".into());
        }
        let status = child.wait().map_err(|e| e.to_string())?;
        if !status.success() {
            return Err(format!("Pixi inspection failed with {status}"));
        }
        String::from_utf8(bytes).map_err(|_| "Pixi inspection output is not UTF-8".into())
    }
    fn control_args(&self, action: &str, name: &str) -> Vec<OsString> {
        [
            OsString::from(action),
            "--no-config".into(),
            "--manifest-path".into(),
            self.workspace().join("pyproject.toml").into_os_string(),
            "--environment".into(),
            name.into(),
        ]
        .into()
    }
    fn verify_workspace(&self) -> Result<(), String> {
        for (name, expected) in [("pyproject.toml", MANIFEST), ("pixi.lock", LOCK)] {
            let path = self.workspace().join(name);
            let meta = fs::symlink_metadata(&path)
                .map_err(|e| format!("workspace is unavailable: {e}; run tools setup"))?;
            if !meta.is_file()
                || meta.len() != expected.len() as u64
                || fs::read(path).map_err(|e| e.to_string())? != expected
            {
                return Err(
                    "bundled workspace contents changed; refusing to use a stale or modified lock"
                        .into(),
                );
            }
        }
        Ok(())
    }
    fn publish_workspace(&self) -> Result<(), String> {
        let workspace = self.workspace();
        if workspace.exists() {
            return self.verify_workspace();
        }
        fs::create_dir_all(workspace.parent().ok_or("missing workspace parent")?)
            .map_err(|e| e.to_string())?;
        let stage = tempfile::tempdir_in(workspace.parent().ok_or("missing workspace parent")?)
            .map_err(|e| e.to_string())?;
        let package = stage.path().join("workspace");
        fs::create_dir(&package).map_err(|e| e.to_string())?;
        for (name, bytes) in [("pyproject.toml", MANIFEST), ("pixi.lock", LOCK)] {
            let mut file = File::create(package.join(name)).map_err(|e| e.to_string())?;
            file.write_all(bytes).map_err(|e| e.to_string())?;
            file.sync_all().map_err(|e| e.to_string())?;
        }
        // The operation lock serializes all Rust publications; an existing
        // nonempty Python-published workspace cannot be replaced by rename.
        if let Err(e) = fs::rename(&package, &workspace) {
            if !workspace.is_dir() {
                return Err(e.to_string());
            }
        }
        self.verify_workspace()
    }
    fn operation_lock(&self, create: bool, exclusive: bool) -> Result<File, String> {
        let directory = self.root.join("native-tool-locks");
        if create {
            fs::create_dir_all(&directory).map_err(|e| e.to_string())?;
        }
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(create)
            .truncate(false)
            .open(directory.join(format!("{}.lock", workspace_identity())))
            .map_err(|e| format!("native installation unavailable: {e}; run tools setup"))?;
        let result = if exclusive {
            FileExt::try_lock_exclusive(&file)
        } else {
            FileExt::try_lock_shared(&file)
        };
        if !result.map_err(|e| format!("cannot lock tool workspace: {e}"))? {
            return Err("tool workspace busy; retry after the active setup/run finishes".into());
        }
        Ok(file)
    }
    fn prefix(&self, name: &str) -> Result<PathBuf, String> {
        let output = self.capture([
            "info".into(),
            "--no-config".into(),
            "--json".into(),
            "--manifest-path".into(),
            self.workspace().join("pyproject.toml").into_os_string(),
        ])?;
        let value: serde_json::Value = serde_json::from_str(&output).map_err(|e| e.to_string())?;
        if value["platform"] != "linux-64" || value["version"] != PIXI_VERSION {
            return Err("Pixi inspection reported unexpected platform/version".into());
        }
        let records = value["environments_info"]
            .as_array()
            .ok_or("Pixi did not report environments_info")?;
        let matching: Vec<_> = records
            .iter()
            .filter(|record| record["name"] == name)
            .collect();
        if matching.len() != 1 {
            return Err("Pixi did not report one selected environment".into());
        }
        let prefix = PathBuf::from(
            matching[0]["prefix"]
                .as_str()
                .ok_or("Pixi did not report a prefix")?,
        );
        // This bridge only operates the conventional bundled workspace. Treat
        // a redirected prefix as unsupported, rather than adopting outside data.
        if prefix != self.workspace().join(".pixi/envs").join(name) {
            return Err("Pixi reported a prefix outside the selected bundled workspace".into());
        }
        if canonical_future_path(&prefix)? != prefix {
            return Err(
                "Pixi prefix resolves through a redirected or symbolic-link location".into(),
            );
        }
        Ok(prefix)
    }
    fn receipt_path(&self, name: &str) -> PathBuf {
        self.workspace().join(format!(".biov-native-{name}.json"))
    }
    fn read_receipt(&self, name: &str) -> Result<Receipt, String> {
        let path = self.receipt_path(name);
        let meta = fs::symlink_metadata(&path).map_err(|_| {
            "no successful native locked setup receipt; run tools setup".to_string()
        })?;
        if !meta.is_file() || meta.len() > MAX_RECORD {
            return Err("native setup receipt is invalid".into());
        }
        let receipt: Receipt = serde_json::from_slice(&fs::read(path).map_err(|e| e.to_string())?)
            .map_err(|e| format!("invalid native setup receipt: {e}"))?;
        let expected =
            self.expected_receipt(name, self.workspace().join(".pixi/envs").join(name))?;
        if receipt != expected {
            return Err("native setup receipt does not match the selected bundled lock/platform; run tools setup".into());
        }
        let marker = receipt.prefix.join("conda-meta/pixi");
        let metadata = fs::symlink_metadata(&marker)
            .map_err(|_| "installed prefix is unavailable; run tools setup".to_string())?;
        if !metadata.is_file() || metadata.len() > MAX_RECORD {
            return Err("Pixi prefix marker is invalid".into());
        }
        let value: serde_json::Value =
            serde_json::from_slice(&fs::read(marker).map_err(|e| e.to_string())?)
                .map_err(|e| e.to_string())?;
        if value["environment_name"] != name
            || value["pixi_version"] != PIXI_VERSION
            || value["resolved_platform"]["subdir"] != "linux-64"
            || value["manifest_path"].as_str() != self.workspace().join("pyproject.toml").to_str()
        {
            return Err("Pixi prefix marker does not match the selected environment".into());
        }
        Ok(receipt)
    }
    fn expected_receipt(&self, name: &str, prefix: PathBuf) -> Result<Receipt, String> {
        let marker = prefix.join("conda-meta/pixi");
        let metadata = fs::symlink_metadata(&marker)
            .map_err(|e| format!("installed Pixi marker is unavailable: {e}"))?;
        if !metadata.is_file() || metadata.len() > MAX_RECORD {
            return Err("installed Pixi marker is invalid".into());
        }
        let prefix_marker_sha256 = digest(&fs::read(marker).map_err(|e| e.to_string())?);
        Ok(Receipt {
            format_version: 1,
            environment: name.into(),
            workspace_sha256: workspace_identity(),
            lock_sha256: digest(LOCK),
            manager_version: PIXI_VERSION.into(),
            platform: "linux-64".into(),
            prefix,
            prefix_marker_sha256,
        })
    }
    fn entrypoint(&self, receipt: &Receipt) -> Result<PathBuf, String> {
        let entrypoint = receipt.prefix.join("bin").join(&receipt.environment);
        let actual = fs::canonicalize(&entrypoint).map_err(|e| format!("selected installed executable is unavailable: {e}; run tools setup after inspecting the prefix"))?;
        let prefix = fs::canonicalize(&receipt.prefix).map_err(|e| e.to_string())?;
        if prefix != receipt.prefix {
            return Err(
                "selected installed prefix resolves outside its recorded bundled location".into(),
            );
        }
        if !actual.starts_with(&prefix) || !actual.is_file() {
            return Err(
                "selected entry point is not an executable file inside the installed prefix".into(),
            );
        }
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            if fs::metadata(&actual)
                .map_err(|e| e.to_string())?
                .permissions()
                .mode()
                & 0o111
                == 0
            {
                return Err("selected installed entry point has no executable permission".into());
            }
        }
        Ok(entrypoint)
    }
    pub fn inspect(&self, name: &str) -> Result<Inspection, String> {
        Self::supported(name)?;
        let checked = self
            .verify_workspace()
            .and_then(|_| self.read_receipt(name))
            .and_then(|receipt| {
                self.entrypoint(&receipt)?;
                Ok(receipt)
            });
        let (status, detail, receipt) = match checked { Ok(receipt) => ("setup_recorded", "Successful locked setup and matching Pixi prefix marker are recorded; package integrity and scientific readiness are not independently verified".into(), Some(receipt)), Err(error) => ("unavailable", error, None) };
        Ok(Inspection {
            environment: name.into(),
            source: "bundled-pixi-lock",
            ownership: "biov-bundled-workspace",
            workspace: self.workspace(),
            workspace_sha256: workspace_identity(),
            lock_sha256: digest(LOCK),
            manager: self.pixi.clone(),
            required_manager_version: PIXI_VERSION,
            platform: "linux-64",
            status,
            detail,
            receipt,
        })
    }
    pub fn setup(&self, name: &str) -> Result<Receipt, String> {
        Self::supported(name)?;
        self.check_manager()?;
        if let Ok(_guard) = self.operation_lock(false, false) {
            if self.verify_workspace().is_ok() {
                if let Ok(receipt) = self.read_receipt(name) {
                    if self.prefix(name)? != receipt.prefix {
                        return Err("current Pixi prefix differs from native setup receipt".into());
                    }
                    self.entrypoint(&receipt)?;
                    return Ok(receipt);
                }
            }
        }
        let _guard = self.operation_lock(true, true)?;
        self.publish_workspace()?;
        if let Ok(receipt) = self.read_receipt(name) {
            self.entrypoint(&receipt)?;
            return Ok(receipt);
        }
        // Invalidate before an in-place native-manager operation: a failed or
        // interrupted install must never retain a successful readiness claim.
        match fs::remove_file(self.receipt_path(name)) {
            Ok(()) => (),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => (),
            Err(e) => return Err(e.to_string()),
        }
        // --no-config still permits project-local Pixi configuration. Reject
        // redirected prefixes before install can write outside this workspace.
        let prefix = self.prefix(name)?;
        let mut args = self.control_args("install", name);
        args.push("--locked".into());
        let status = Command::new(&self.pixi)
            .args(args)
            .stdin(Stdio::inherit())
            .stdout(Stdio::inherit())
            .stderr(Stdio::inherit())
            .status()
            .map_err(|e| e.to_string())?;
        if !status.success() {
            return Err(format!("locked Pixi setup failed with {status}; native readiness is unavailable; existing files were not removed by BioV"));
        }
        self.verify_workspace()?;
        if self.prefix(name)? != prefix {
            return Err("Pixi prefix changed during setup".into());
        }
        let receipt = self.expected_receipt(name, prefix)?;
        let mut staged =
            tempfile::NamedTempFile::new_in(self.workspace()).map_err(|e| e.to_string())?;
        serde_json::to_writer_pretty(&mut staged, &receipt).map_err(|e| e.to_string())?;
        staged.write_all(b"\n").map_err(|e| e.to_string())?;
        staged.as_file().sync_all().map_err(|e| e.to_string())?;
        staged
            .persist(self.receipt_path(name))
            .map_err(|e| e.to_string())?;
        let receipt = self.read_receipt(name)?;
        self.entrypoint(&receipt)?;
        Ok(receipt)
    }
    /// Pixi activates the environment, then runs its exact prefix executable.
    /// No named task interpolation or selected-entry host-PATH fallback occurs.
    pub fn execute(
        &self,
        name: &str,
        arguments: &[OsString],
        cwd: Option<&Path>,
    ) -> Result<ExitStatus, String> {
        Self::supported(name)?;
        if arguments.iter().any(|argument| argument.to_str().is_none()) {
            return Err("Pixi native arguments must be UTF-8".into());
        }
        self.check_manager()?;
        let _guard = self.operation_lock(false, false)?;
        self.verify_workspace()?;
        let receipt = self.read_receipt(name)?;
        if self.prefix(name)? != receipt.prefix {
            return Err("current Pixi prefix differs from native setup receipt".into());
        }
        let entrypoint = self.entrypoint(&receipt)?;
        let mut args = self.control_args("run", name);
        args.push("--as-is".into());
        args.push("--executable".into());
        args.push("--".into());
        args.push(entrypoint.into_os_string());
        args.extend_from_slice(arguments);
        let mut command = Command::new(&self.pixi);
        command
            .args(args)
            .stdin(Stdio::inherit())
            .stdout(Stdio::inherit())
            .stderr(Stdio::inherit());
        if let Some(cwd) = cwd {
            if !cwd.is_dir() {
                return Err("execution cwd must be an existing directory".into());
            }
            command.current_dir(cwd);
        }
        command
            .status()
            .map_err(|e| format!("Pixi execution failed: {e}"))
    }
}

// Pixi canonicalizes the workspace before reporting its prefix. Resolve the
// existing ancestor now without creating anything during inspection.
fn canonical_future_path(path: &Path) -> Result<PathBuf, String> {
    match fs::canonicalize(path) {
        Ok(path) => Ok(path),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            let name = path.file_name().ok_or("cannot resolve tool root")?;
            let parent = path.parent().ok_or("cannot resolve tool root parent")?;
            Ok(canonical_future_path(parent)?.join(name))
        }
        Err(error) => Err(format!("cannot resolve tool root: {error}")),
    }
}
