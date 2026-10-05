//! Local locked-tool bridge. Pixi owns resolution, installation and activation.
//! This first migration slice supports two preparation-free bundled tools.
use fs4::fs_std::FileExt;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    ffi::OsString,
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    process::{Command, ExitStatus, Stdio},
};

mod installed;
pub mod model_bundle;
pub mod model_resources;

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
        self.root.join("workspaces").join(self.identity())
    }
    fn identity(&self) -> String {
        workspace_identity()
    }
    fn lock_digest(&self) -> String {
        digest(LOCK)
    }
    fn supported(name: &str) -> Result<(), String> {
        if SUPPORTED_TOOLS.contains(&name) {
            Ok(())
        } else {
            Err(format!("unsupported native tool {name:?}; initial tools are samtools and goatools; use biov python for other environments"))
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
        self.capture_in(args, None)
    }
    fn capture_in(
        &self,
        args: impl IntoIterator<Item = OsString>,
        cwd: Option<&Path>,
    ) -> Result<String, String> {
        let mut command = Command::new(&self.pixi);
        command
            .args(args)
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit());
        if let Some(cwd) = cwd {
            command.current_dir(cwd);
        }
        let mut child = command.spawn().map_err(|e| {
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
            return Err(format!("Pixi query failed with {status}"));
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
            let meta = fs::symlink_metadata(&path).map_err(|e| {
                format!("workspace is unavailable: {e}; run biov tools exec <tool>")
            })?;
            if !meta.is_file()
                || meta.len() > 16 * 1024 * 1024
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
            .open(directory.join(format!("{}.lock", self.identity())))
            .map_err(|e| {
                format!("native installation unavailable: {e}; run biov tools exec <tool>")
            })?;
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
            "no successful native locked setup receipt; run biov tools exec <tool>".to_string()
        })?;
        if !meta.is_file() || meta.len() > MAX_RECORD {
            return Err("native setup receipt is invalid".into());
        }
        let receipt: Receipt = serde_json::from_slice(&fs::read(path).map_err(|e| e.to_string())?)
            .map_err(|e| format!("invalid native setup receipt: {e}"))?;
        let expected =
            self.expected_receipt(name, self.workspace().join(".pixi/envs").join(name))?;
        if receipt != expected {
            return Err("native setup receipt does not match the selected bundled lock/platform; run biov tools exec <tool>".into());
        }
        let marker = receipt.prefix.join("conda-meta/pixi");
        let metadata = fs::symlink_metadata(&marker).map_err(|_| {
            "installed prefix is unavailable; run biov tools exec <tool>".to_string()
        })?;
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
            workspace_sha256: self.identity(),
            lock_sha256: self.lock_digest(),
            manager_version: PIXI_VERSION.into(),
            platform: "linux-64".into(),
            prefix,
            prefix_marker_sha256,
        })
    }
    fn entrypoint_parent(&self, prefix: &Path) -> Result<(), String> {
        if canonical_future_path(prefix)? != prefix {
            return Err("selected installed prefix is redirected".into());
        }
        let bin = prefix.join("bin");
        if canonical_future_path(&bin)? != bin {
            return Err("selected executable bin directory is redirected".into());
        }
        match fs::metadata(&bin) {
            Ok(metadata) if !metadata.is_dir() => {
                Err("selected executable bin path is not a directory".into())
            }
            Ok(_) => Ok(()),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
            Err(error) => Err(format!(
                "cannot inspect selected executable bin directory: {error}"
            )),
        }
    }
    fn entrypoint_location(&self, prefix: &Path, name: &str) -> Result<(), String> {
        self.entrypoint_parent(prefix)?;
        let entrypoint = prefix.join("bin").join(name);
        match fs::symlink_metadata(&entrypoint) {
            Ok(metadata) if metadata.file_type().is_symlink() => {
                let target = fs::canonicalize(&entrypoint).map_err(|e| {
                    format!("selected executable symlink is unavailable or unsafe: {e}")
                })?;
                if !target.starts_with(prefix) {
                    return Err(
                        "selected executable symlink resolves outside the installed prefix".into(),
                    );
                }
                Ok(())
            }
            Ok(_) => Ok(()),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
            Err(error) => Err(format!(
                "cannot inspect selected executable location: {error}"
            )),
        }
    }
    fn entrypoint(&self, receipt: &Receipt) -> Result<PathBuf, String> {
        self.entrypoint_location(&receipt.prefix, &receipt.environment)?;
        let entrypoint = receipt.prefix.join("bin").join(&receipt.environment);
        let actual = fs::canonicalize(&entrypoint).map_err(|e| format!("selected installed executable is unavailable: {e}; run biov tools exec <tool> after inspecting the prefix"))?;
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
            workspace_sha256: self.identity(),
            lock_sha256: self.lock_digest(),
            manager: self.pixi.clone(),
            required_manager_version: PIXI_VERSION,
            platform: "linux-64",
            status,
            detail,
            receipt,
        })
    }
    fn missing_entrypoint(&self, prefix: &Path, name: &str) -> Result<bool, String> {
        self.entrypoint_location(prefix, name)?;
        match fs::symlink_metadata(prefix.join("bin").join(name)) {
            Ok(_) => Ok(false),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(true),
            Err(error) => Err(format!(
                "cannot inspect selected installed executable: {error}"
            )),
        }
    }
    fn reusable_receipt(&self, name: &str) -> Result<Option<Receipt>, String> {
        let receipt = match self.read_receipt(name) {
            Ok(receipt) => receipt,
            Err(_) => return Ok(None),
        };
        if self.prefix(name)? != receipt.prefix {
            return Err("current Pixi prefix differs from native setup receipt".into());
        }
        match self.entrypoint(&receipt) {
            Ok(_) => Ok(Some(receipt)),
            Err(_) if self.missing_entrypoint(&receipt.prefix, name)? => Ok(None),
            Err(error) => Err(error),
        }
    }
    fn install(&self, action: &str, name: &str) -> Result<(), String> {
        let mut args = self.control_args(action, name);
        args.push("--locked".into());
        let status = Command::new(&self.pixi)
            .args(args)
            .stdin(Stdio::inherit())
            .stdout(Stdio::inherit())
            .stderr(Stdio::inherit())
            .status()
            .map_err(|e| e.to_string())?;
        if !status.success() {
            return Err(format!("locked Pixi {action} failed with {status}; native readiness is unavailable; existing files were not removed by BioV"));
        }
        self.verify_workspace()
    }
    pub fn setup(&self, name: &str) -> Result<Receipt, String> {
        Self::supported(name)?;
        self.check_manager()?;
        if let Ok(_guard) = self.operation_lock(false, false) {
            if self.verify_workspace().is_ok() {
                if let Some(receipt) = self.reusable_receipt(name)? {
                    return Ok(receipt);
                }
            }
        }
        let _guard = self.operation_lock(true, true)?;
        self.publish_workspace()?;
        // Apply identical reuse/routing checks after obtaining the exclusive
        // lock: another setup could have completed since the shared fast path.
        if let Some(receipt) = self.reusable_receipt(name)? {
            return Ok(receipt);
        }
        let prefix = self.prefix(name)?;
        self.entrypoint_location(&prefix, name)?;
        match fs::remove_file(self.receipt_path(name)) {
            Ok(()) => (),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => (),
            Err(e) => return Err(e.to_string()),
        }
        self.install("install", name)?;
        if self.prefix(name)? != prefix {
            return Err("Pixi prefix changed during setup".into());
        }
        // A normal manager install can regard an existing package record as
        // current even after its executable was removed. Only a demonstrably
        // absent entry point triggers explicit locked manager reinstallation.
        if self.missing_entrypoint(&prefix, name)? {
            self.install("reinstall", name)?;
            if self.prefix(name)? != prefix {
                return Err("Pixi prefix changed during repair".into());
            }
        }
        let receipt = self.expected_receipt(name, prefix)?;
        self.entrypoint(&receipt)?;
        let mut staged =
            tempfile::NamedTempFile::new_in(self.workspace()).map_err(|e| e.to_string())?;
        serde_json::to_writer_pretty(&mut staged, &receipt).map_err(|e| e.to_string())?;
        staged.write_all(b"\n").map_err(|e| e.to_string())?;
        staged.as_file().sync_all().map_err(|e| e.to_string())?;
        staged
            .persist(self.receipt_path(name))
            .map_err(|e| e.to_string())?;
        self.read_receipt(name)
    }
    fn activation(
        &self,
        name: &str,
        receipt: &Receipt,
        cwd: Option<&Path>,
    ) -> Result<BTreeMap<String, String>, String> {
        let mut args = self.control_args("shell-hook", name);
        args.extend(["--as-is".into(), "--json".into()]);
        let output = self.capture_in(args, cwd)?;
        // Pixi 0.81.0 already evaluates manifest/conda activation scripts while
        // computing these variables. Its activation_scripts field is informational,
        // not pending work: sourcing it again would repeat script side effects.
        #[derive(Deserialize)]
        struct Activation {
            environment_variables: BTreeMap<String, String>,
        }
        let activation: Activation = serde_json::from_str(&output)
            .map_err(|e| format!("invalid Pixi activation JSON: {e}"))?;
        let variables = activation.environment_variables;
        for (key, value) in &variables {
            if key.is_empty() || key.contains(['=', '\0']) || value.contains('\0') {
                return Err("Pixi activation contains an invalid environment variable".into());
            }
        }
        if variables.get("CONDA_PREFIX").map(String::as_str) != receipt.prefix.to_str()
            || variables.get("PIXI_ENVIRONMENT_NAME").map(String::as_str) != Some(name)
            || variables.get("PIXI_PROJECT_MANIFEST").map(String::as_str)
                != self.workspace().join("pyproject.toml").to_str()
            || variables.get("PIXI_PROJECT_ROOT").map(String::as_str) != self.workspace().to_str()
        {
            return Err(
                "Pixi activation does not match the selected bundled prefix/environment/manifest"
                    .into(),
            );
        }
        let path = variables
            .get("PATH")
            .ok_or("Pixi activation did not report PATH")?;
        if std::env::split_paths(path).next() != Some(receipt.prefix.join("bin")) {
            return Err("Pixi activation PATH does not start with the selected prefix/bin".into());
        }
        Ok(variables)
    }
    /// Apply the manager's typed JSON activation, then launch the exact prefix
    /// executable with OS argv. No Pixi task/command-string parsing occurs.
    pub fn execute(
        &self,
        name: &str,
        arguments: &[OsString],
        cwd: Option<&Path>,
    ) -> Result<ExitStatus, String> {
        self.execute_using(name, arguments, cwd, |command| {
            command
                .status()
                .map_err(|e| format!("native executable launch failed: {e}"))
        })
    }
    /// CLI-owned cooperative interrupt policy. This explicitly installs Unix
    /// process-wide signal handlers, forwarding INT/TERM/HUP during execution and
    /// using default dispositions afterward. Embedders should use `execute` and
    /// retain their own process signal policy instead.
    pub fn execute_with_interrupt_forwarding(
        &self,
        name: &str,
        arguments: &[OsString],
        cwd: Option<&Path>,
    ) -> Result<ExitStatus, String> {
        self.execute_using(name, arguments, cwd, crate::execution::run)
    }
    fn execute_using(
        &self,
        name: &str,
        arguments: &[OsString],
        cwd: Option<&Path>,
        run: impl FnOnce(&mut Command) -> Result<ExitStatus, String>,
    ) -> Result<ExitStatus, String> {
        Self::supported(name)?;
        if arguments.iter().any(|argument| argument.to_str().is_none()) {
            return Err("native tool arguments must be UTF-8".into());
        }
        self.check_manager()?;
        let _guard = self.operation_lock(false, false)?;
        self.verify_workspace()?;
        let receipt = self.read_receipt(name)?;
        if self.prefix(name)? != receipt.prefix {
            return Err("current Pixi prefix differs from native setup receipt".into());
        }
        let entrypoint = self.entrypoint(&receipt)?;
        if let Some(cwd) = cwd {
            if !cwd.is_dir() {
                return Err("execution cwd must be an existing directory".into());
            }
        }
        let activation = self.activation(name, &receipt, cwd)?;
        let mut command = Command::new(entrypoint);
        command
            .args(arguments)
            .envs(activation)
            .stdin(Stdio::inherit())
            .stdout(Stdio::inherit())
            .stderr(Stdio::inherit());
        if let Some(cwd) = cwd {
            command.current_dir(cwd);
        }
        run(&mut command)
    }
}

mod execution;

// Pixi canonicalizes the workspace before reporting its prefix. Resolve the
// existing ancestor now without creating anything during inspection.
fn canonical_future_path(path: &Path) -> Result<PathBuf, String> {
    match fs::canonicalize(path) {
        Ok(path) => Ok(path),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            if path.components().next_back() == Some(std::path::Component::ParentDir) {
                let parent =
                    canonical_future_path(path.parent().ok_or("cannot resolve parent component")?)?;
                return parent
                    .parent()
                    .map(Path::to_path_buf)
                    .ok_or_else(|| "path escapes filesystem root".into());
            }
            let name = path.file_name().ok_or("cannot resolve tool root")?;
            let parent = path.parent().ok_or("cannot resolve tool root parent")?;
            Ok(canonical_future_path(parent)?.join(name))
        }
        Err(error) => Err(format!("cannot resolve tool root: {error}")),
    }
}

#[cfg(all(test, target_os = "linux", target_arch = "x86_64"))]
mod real_pixi_tests;
