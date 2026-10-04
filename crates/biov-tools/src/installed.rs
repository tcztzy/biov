//! Convenience installations belong to upstream global/tool backends. The
//! bundled locked workspace remains separate for scientific `tools exec`.
use super::*;

const SAMTOOLS: &str = "samtools==1.24";
const GOATOOLS: &str = "goatools==1.6.5";
const STATSMODELS: &str = "statsmodels==0.14.6";
// The pinned GOATOOLS wheel's standard console_scripts entry points. uv owns
// their creation and removal; this list only bounds collision preflight.
const GO_COMMANDS: &[&str] = &[
    "goatools",
    "find_enrichment.py",
    "go_plot.py",
    "map_to_slim.py",
    "ncbi_gene_results_to_python.py",
    "plot_go_term.py",
    "wr_hier.py",
];

fn direct_location(path: &Path) -> Result<(), String> {
    if fs::symlink_metadata(path).is_ok_and(|m| m.file_type().is_symlink())
        || canonical_future_path(path)? != path
    {
        return Err(format!(
            "backend location is redirected: {}",
            path.display()
        ));
    }
    Ok(())
}

fn exists(path: &Path) -> Result<bool, String> {
    match fs::symlink_metadata(path) {
        Ok(_) => Ok(true),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(false),
        Err(e) => Err(format!("cannot inspect {}: {e}", path.display())),
    }
}

fn run_backend(command: &mut Command, label: &str) -> Result<(), String> {
    let status = command
        .stdin(Stdio::inherit())
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .status()
        .map_err(|e| format!("cannot run {label}: {e}"))?;
    if !status.success() {
        return Err(format!(
            "{label} failed with {status}; inspect its native state before retrying"
        ));
    }
    Ok(())
}

fn capture_backend(command: &mut Command, label: &str) -> Result<String, String> {
    let mut child = command
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .map_err(|e| format!("cannot run {label}: {e}"))?;
    let mut bytes = Vec::new();
    let read = child
        .stdout
        .take()
        .ok_or("missing backend stdout")?
        .take(1024 * 1024 + 1)
        .read_to_end(&mut bytes);
    if read.is_err() || bytes.len() > 1024 * 1024 {
        let _ = child.kill();
        let _ = child.wait();
        return Err(format!("{label} output exceeds 1 MiB or could not be read"));
    }
    let status = child.wait().map_err(|e| e.to_string())?;
    if !status.success() {
        return Err(format!("{label} failed with {status}"));
    }
    String::from_utf8(bytes).map_err(|_| format!("{label} output is not UTF-8"))
}

impl ToolStore {
    fn global_pixi_home(&self) -> PathBuf {
        self.root.join("pixi-global")
    }
    fn global_pixi_manifest(&self) -> PathBuf {
        self.global_pixi_home().join("manifests/pixi-global.toml")
    }
    fn global_uv_home(&self) -> PathBuf {
        self.root.join("uv-tools")
    }
    /// Upstream-owned command directories. BioV never edits shell profiles.
    pub fn bin_paths(&self) -> Vec<PathBuf> {
        vec![
            self.global_pixi_home().join("bin"),
            self.global_uv_home().join("bin"),
        ]
    }
    fn pixi_global_command(&self, action: &str) -> Result<Command, String> {
        let home = self.global_pixi_home();
        for path in [
            home.clone(),
            home.join("bin"),
            home.join("envs"),
            home.join("manifests"),
            home.join("cache"),
        ] {
            direct_location(&path)?;
        }
        let mut command = Command::new(&self.pixi);
        command
            .args(["global", action, "--color", "never"])
            .env("PIXI_HOME", &home)
            .env("PIXI_CACHE_DIR", home.join("cache"))
            .env("XDG_CONFIG_HOME", home.join("config"))
            .env_remove("PIXI_CONFIG_FILE");
        Ok(command)
    }
    fn uv_command(&self, uv: Option<&Path>) -> Result<Command, String> {
        let uv = uv
            .map(Path::to_path_buf)
            .or_else(|| std::env::var_os("BIOV_UV_BIN").map(PathBuf::from))
            .unwrap_or_else(|| PathBuf::from("uv"));
        let home = self.global_uv_home();
        for path in [
            home.clone(),
            home.join("tools"),
            home.join("bin"),
            home.join("cache"),
            home.join("python"),
        ] {
            direct_location(&path)?;
        }
        let mut command = Command::new(uv);
        command
            .env("UV_TOOL_DIR", home.join("tools"))
            .env("UV_TOOL_BIN_DIR", home.join("bin"))
            .env("UV_CACHE_DIR", home.join("cache"))
            .env("UV_PYTHON_INSTALL_DIR", home.join("python"))
            .env_remove("UV_CONFIG_FILE");
        Ok(command)
    }
    fn require_pixi_manifest(&self) -> Result<(), String> {
        let path = self.global_pixi_manifest();
        direct_location(&path)?;
        let meta = fs::symlink_metadata(&path)
            .map_err(|e| format!("isolated Pixi global manifest unavailable: {e}"))?;
        if !meta.is_file() {
            return Err("isolated Pixi global manifest must be a regular file".into());
        }
        Ok(())
    }
    fn prepare_pixi_manifest(&self) -> Result<(), String> {
        let manifest = self.global_pixi_manifest();
        if exists(&manifest)? {
            return self.require_pixi_manifest();
        }
        let parent = manifest.parent().ok_or("missing Pixi manifest parent")?;
        direct_location(parent)?;
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        // Pixi 0.81 searches both PIXI_HOME and XDG config for a pre-existing
        // manifest. Publish its standard native TOML here before invoking it,
        // so it cannot adopt unrelated user-global installations.
        let mut stage = tempfile::NamedTempFile::new_in(parent).map_err(|e| e.to_string())?;
        stage
            .write_all(b"version = 1\n")
            .map_err(|e| e.to_string())?;
        stage.as_file().sync_all().map_err(|e| e.to_string())?;
        match stage.persist_noclobber(&manifest) {
            Ok(_) => (),
            Err(_) if exists(&manifest)? => (),
            Err(e) => return Err(e.to_string()),
        }
        self.require_pixi_manifest()
    }
    fn pixi_samtools_selected(&self) -> Result<bool, String> {
        if !exists(&self.global_pixi_manifest())? {
            return Ok(false);
        }
        self.require_pixi_manifest()?;
        let mut command = self.pixi_global_command("list")?;
        command.args(["--json", "--offline"]);
        // This is Pixi's documented CLI output, not a private receipt schema.
        let output = capture_backend(&mut command, "Pixi global inventory")?;
        let environments: Vec<serde_json::Value> = serde_json::from_str(&output)
            .map_err(|e| format!("invalid Pixi global inventory: {e}"))?;
        let selected: Vec<_> = environments
            .iter()
            .filter(|r| r["name"] == "samtools")
            .collect();
        if selected.is_empty() {
            return Ok(false);
        }
        if selected.len() != 1
            || selected[0]["dependencies"]
                != serde_json::json!([{"name":"samtools","version":"1.24"}])
            || selected[0]["exposed"]
                != serde_json::json!([{"exposed_name":"samtools","executable":"samtools"}])
        {
            return Err("isolated Samtools environment has an unexpected upstream selection; inspect Pixi global state".into());
        }
        Ok(true)
    }
    fn preflight_pixi_command(&self) -> Result<(), String> {
        let bin = self.global_pixi_home().join("bin");
        let target = bin.join("samtools");
        if !exists(&target)? {
            return Ok(());
        }
        let shared = bin.join("trampoline_configuration/trampoline_bin");
        direct_location(&shared)?;
        for path in [&target, &shared] {
            let meta = fs::symlink_metadata(path).map_err(|e| e.to_string())?;
            if !meta.is_file() || meta.len() > 16 * 1024 * 1024 {
                return Err(format!(
                    "refusing modified or foreign command: {}",
                    target.display()
                ));
            }
        }
        if !self.pixi_samtools_selected()?
            || fs::read(&target).map_err(|e| e.to_string())?
                != fs::read(&shared).map_err(|e| e.to_string())?
        {
            return Err(format!(
                "refusing modified or foreign command: {}",
                target.display()
            ));
        }
        Ok(())
    }
    fn preflight_uv_commands(&self) -> Result<(), String> {
        let home = self.global_uv_home();
        let prefix = home.join("tools/goatools");
        direct_location(&prefix)?;
        direct_location(&prefix.join("bin"))?;
        for name in GO_COMMANDS {
            let target = home.join("bin").join(name);
            if !exists(&target)? {
                continue;
            }
            let expected = prefix.join("bin").join(name);
            let meta = fs::symlink_metadata(&target).map_err(|e| e.to_string())?;
            if !meta.file_type().is_symlink()
                || fs::canonicalize(&target).ok() != Some(expected.clone())
                || !expected.is_file()
            {
                return Err(format!(
                    "refusing modified or foreign command: {}",
                    target.display()
                ));
            }
        }
        Ok(())
    }
    /// Install a directly usable upstream command. Direct package versions are
    /// pinned; this convenience route has no frozen transitive scientific lock.
    pub fn install_global(
        &self,
        name: &str,
        python: Option<&Path>,
        uv: Option<&Path>,
    ) -> Result<(), String> {
        Self::supported(name)?;
        if name == "samtools" {
            self.check_manager()?;
            let mut command = self.pixi_global_command("install")?;
            self.preflight_pixi_command()?;
            self.prepare_pixi_manifest()?;
            command.args([
                SAMTOOLS,
                "--channel",
                "conda-forge",
                "--channel",
                "bioconda",
                "--platform",
                "linux-64",
                "--environment",
                "samtools",
                "--expose",
                "samtools",
                "--no-shortcuts",
            ]);
            run_backend(&mut command, "Pixi global Samtools install")
        } else {
            let python = python.ok_or(
                "GOATOOLS installation requires the Python interpreter paired with installed BioV",
            )?;
            if !python.is_file() {
                return Err("paired Python interpreter is unavailable".into());
            }
            let mut command = self.uv_command(uv)?;
            self.preflight_uv_commands()?;
            command
                .args([
                    "tool",
                    "install",
                    GOATOOLS,
                    "--with",
                    STATSMODELS,
                    "--python",
                ])
                .arg(python)
                .args(["--no-python-downloads", "--no-config"]);
            run_backend(&mut command, "uv tool GOATOOLS install")
        }
    }
    /// Display upstream inventories only from the dedicated backend roots.
    /// No private uv receipt parsing or normalized JSON contract is provided.
    pub fn list_global(&self, uv: Option<&Path>) -> Result<String, String> {
        let mut sections = Vec::new();
        if exists(&self.global_pixi_manifest())? {
            self.check_manager()?;
            self.require_pixi_manifest()?;
            let mut command = self.pixi_global_command("list")?;
            command.arg("--offline");
            sections.push(format!(
                "Pixi global ({})\n{}",
                self.global_pixi_home().display(),
                capture_backend(&mut command, "Pixi global list")?
            ));
        }
        let tools = self.global_uv_home().join("tools");
        if exists(&tools)? {
            let mut command = self.uv_command(uv)?;
            command.args([
                "tool",
                "list",
                "--show-paths",
                "--show-with",
                "--show-version-specifiers",
                "--offline",
                "--no-config",
                "--color",
                "never",
            ]);
            sections.push(format!(
                "uv tool ({})\n{}",
                tools.display(),
                capture_backend(&mut command, "uv tool list")?
            ));
        }
        if sections.is_empty() {
            return Ok("No native tools installed in this environment root.\n".into());
        }
        Ok(sections.join("\n"))
    }
    /// Delegate removal of the named backend environment and commands. Locked
    /// scientific workspaces, caches, biological data and outputs are untouched.
    pub fn uninstall_global(&self, name: &str, uv: Option<&Path>) -> Result<(), String> {
        Self::supported(name)?;
        if name == "samtools" {
            self.check_manager()?;
            self.require_pixi_manifest()?;
            let mut command = self.pixi_global_command("uninstall")?;
            self.preflight_pixi_command()?;
            if !self.pixi_samtools_selected()? {
                return Err("Samtools is not installed in the isolated Pixi global root".into());
            }
            command.arg("samtools");
            run_backend(&mut command, "Pixi global Samtools uninstall")
        } else {
            let tools = self.global_uv_home().join("tools");
            if !exists(&tools)? {
                return Err("GOATOOLS is not installed in the isolated uv tool root".into());
            }
            let mut command = self.uv_command(uv)?;
            self.preflight_uv_commands()?;
            command.args([
                "tool",
                "uninstall",
                "goatools",
                "--no-config",
                "--no-python-downloads",
            ]);
            run_backend(&mut command, "uv tool GOATOOLS uninstall")
        }
    }
}
