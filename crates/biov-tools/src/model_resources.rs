//! Model files remain data. Official hf owns Hub access; uv owns fallback tools.
use crate::model_bundle::{inspect_bundle, validate_selection, write_bundle, ModelResourceRecord};
use serde::Serialize;
use std::{
    fs,
    io::Read,
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

pub const HF_PACKAGE: &str = "huggingface-hub==2.1.1";
const HF_FALLBACK_VERSION: &str = "2.1.1";
const MAX_QUERY: u64 = 128 * 1024;

pub struct ModelDownload {
    pub repository: String,
    pub revision: String,
    pub files: Vec<String>,
    pub local_dir: PathBuf,
    pub hf: Option<PathBuf>,
    pub uv: Option<PathBuf>,
    pub no_install: bool,
}

#[derive(Serialize)]
pub struct ModelResourceResult {
    pub status: &'static str,
    pub local_dir: PathBuf,
    pub resource: ModelResourceRecord,
}

enum Client {
    Installed(PathBuf),
    Uv(PathBuf),
}

impl Client {
    fn command(&self) -> Command {
        match self {
            Self::Installed(path) => Command::new(path),
            Self::Uv(path) => {
                let mut command = Command::new(path);
                command.args([
                    "tool",
                    "run",
                    "--no-config",
                    "--no-python-downloads",
                    "--from",
                    HF_PACKAGE,
                    "--with",
                    "httpx2[socks]",
                    "hf",
                ]);
                command
            }
        }
    }
    fn version(&self) -> Result<String, String> {
        let output = capture(self.command().arg("version"))?;
        let version = parse_version(&output)?;
        if matches!(self, Self::Uv(_)) && version != HF_FALLBACK_VERSION {
            return Err("uv fallback did not supply the requested pinned hf version".into());
        }
        let help = capture(self.command().args(["download", "--help"]))?;
        if ["--revision", "--repo-type", "--local-dir"]
            .iter()
            .any(|flag| !help.contains(flag))
        {
            return Err("hf download does not expose the required upstream options".into());
        }
        Ok(version)
    }
}

fn capture(command: &mut Command) -> Result<String, String> {
    let mut child = command
        .env("NO_COLOR", "1")
        .env("TERM", "dumb")
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .map_err(|e| format!("cannot query upstream hf: {e}"))?;
    let mut bytes = Vec::new();
    let read = child
        .stdout
        .take()
        .ok_or("hf query has no stdout")?
        .take(MAX_QUERY + 1)
        .read_to_end(&mut bytes);
    if read.is_err() || bytes.len() as u64 > MAX_QUERY {
        let _ = child.kill();
        let _ = child.wait();
        return Err("upstream hf query exceeds 128 KiB or cannot be read".into());
    }
    let status = child
        .wait()
        .map_err(|e| format!("cannot reap hf query: {e}"))?;
    if !status.success() {
        return Err(format!("upstream hf query failed with {status}"));
    }
    String::from_utf8(bytes).map_err(|_| "hf query is not UTF-8".into())
}

fn parse_version(output: &str) -> Result<String, String> {
    let mut found = None;
    for line in output.lines().map(str::trim) {
        let value = ["huggingface_hub version:", "version:", "version="]
            .iter()
            .find_map(|prefix| line.strip_prefix(prefix))
            .unwrap_or(line)
            .trim();
        let parts: Vec<_> = value.split('.').collect();
        if parts.len() != 3
            || parts
                .iter()
                .any(|part| part.is_empty() || !part.bytes().all(|byte| byte.is_ascii_digit()))
        {
            continue;
        }
        let major: u64 = parts[0].parse().map_err(|_| "invalid hf version")?;
        let minor: u64 = parts[1].parse().map_err(|_| "invalid hf version")?;
        if major >= 3 || (major == 0 && minor < 34) {
            return Err(
                "hf must be a stable version >=0.34 and <3 with the required download options"
                    .into(),
            );
        }
        if found.replace(value.to_owned()).is_some() {
            return Err("hf version query contains multiple version fields".into());
        }
    }
    found.ok_or_else(|| "unrecognized hf version output; a stable upstream hf is required".into())
}

fn select_client(request: &ModelDownload) -> Result<(Client, String), String> {
    let explicit = request
        .hf
        .clone()
        .or_else(|| std::env::var_os("BIOV_HF_BIN").map(PathBuf::from));
    let candidate = resolve_program(explicit.as_deref().unwrap_or_else(|| Path::new("hf")))
        .map(Client::Installed);
    match candidate.and_then(|candidate| candidate.version().map(|version| (candidate, version))) {
        Ok(client) => return Ok(client),
        Err(error) if explicit.is_some() || request.no_install => return Err(error),
        Err(_) => eprintln!(
            "biov: no compatible hf found on PATH; using uv's pinned on-demand {HF_PACKAGE}"
        ),
    }
    let uv = request
        .uv
        .clone()
        .or_else(|| std::env::var_os("BIOV_UV_BIN").map(PathBuf::from))
        .unwrap_or_else(|| "uv".into());
    let client = Client::Uv(resolve_program(&uv)?);
    let version = client.version()?;
    Ok((client, version))
}

// Resolve once in the caller's cwd. Relative PATH entries and explicit relative
// executable paths retain their original meaning alongside absolute local-dir.
fn resolve_program(program: &Path) -> Result<PathBuf, String> {
    if program.is_absolute() || program.components().count() > 1 {
        return fs::canonicalize(program).map_err(|e| {
            format!(
                "upstream executable {} is unavailable: {e}",
                program.display()
            )
        });
    }
    for directory in std::env::split_paths(&std::env::var_os("PATH").unwrap_or_default()) {
        let candidate = directory.join(program);
        if !candidate.is_file() {
            continue;
        }
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            if fs::metadata(&candidate)
                .map_err(|e| e.to_string())?
                .permissions()
                .mode()
                & 0o111
                == 0
            {
                continue;
            }
        }
        return fs::canonicalize(candidate).map_err(|e| e.to_string());
    }
    Err(format!(
        "upstream executable {} is not on PATH",
        program.display()
    ))
}

fn same_selection(record: &ModelResourceRecord, request: &ModelDownload) -> bool {
    let mut files = request.files.clone();
    files.sort();
    record.repository == request.repository
        && record.revision == request.revision
        && record.selected_files == files
}

pub fn inspect(local_dir: &Path) -> Result<ModelResourceResult, String> {
    let resource = inspect_bundle(local_dir)?;
    let local_dir = fs::canonicalize(local_dir).map_err(|e| e.to_string())?;
    Ok(ModelResourceResult {
        status: "verified",
        local_dir,
        resource,
    })
}

/// Publish a new selected-file bundle, or verify and reuse the exact existing one.
/// No existing destination is adopted, repaired, overwritten or refreshed.
pub fn download(request: &ModelDownload) -> Result<ModelResourceResult, String> {
    validate_selection(&request.repository, &request.revision, &request.files)?;
    if !cfg!(all(target_os = "linux", target_arch = "x86_64")) {
        return Err("model-resource downloads currently support Linux x86_64".into());
    }
    if request.local_dir.to_str().is_none() {
        return Err("model-resource destination must be UTF-8".into());
    }
    let local_dir = if request.local_dir.is_absolute() {
        request.local_dir.clone()
    } else {
        std::env::current_dir()
            .map_err(|e| format!("cannot resolve destination cwd: {e}"))?
            .join(&request.local_dir)
    };
    validate_destination_components(&local_dir)?;
    match fs::symlink_metadata(&local_dir) {
        Ok(metadata) => {
            if !metadata.is_dir() || metadata.file_type().is_symlink() {
                return Err("model-resource destination must be a direct directory".into());
            }
            let mut result = inspect(&local_dir)?;
            if !same_selection(&result.resource, request) {
                return Err("existing model resource differs from the requested repository/revision/files; choose a new directory".into());
            }
            result.status = "reused";
            return Ok(result);
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => (),
        Err(error) => {
            return Err(format!(
                "cannot inspect model-resource destination: {error}"
            ))
        }
    }
    let destination = crate::canonical_future_path(&local_dir)?;
    let parent = destination
        .parent()
        .ok_or("model-resource destination has no parent")?;
    if destination.file_name().is_none() {
        return Err("model-resource destination needs a directory name".into());
    }
    if std::env::var_os("HF_ENDPOINT").is_some_and(|endpoint| {
        endpoint
            .to_str()
            .is_none_or(|text| text.trim_end_matches('/') != "https://huggingface.co")
    }) {
        return Err(
            "custom HF_ENDPOINT is unsupported by this Hugging Face source contract".into(),
        );
    }
    if std::env::var_os("HUGGINGFACE_CO_STAGING").is_some_and(|staging| {
        staging.to_str().is_none_or(|value| {
            // Match upstream Python str.upper(), including Unicode such as yeſ.
            ["1", "ON", "YES", "TRUE"].contains(&value.to_uppercase().as_str())
        })
    }) {
        return Err(
            "HUGGINGFACE_CO_STAGING is unsupported by this public Hugging Face source contract"
                .into(),
        );
    }
    let (client, version) = select_client(request)?;
    fs::create_dir_all(parent).map_err(|e| format!("cannot prepare resource parent: {e}"))?;
    let staging = tempfile::Builder::new()
        .prefix(".biov-model-download-")
        .tempdir_in(parent)
        .map_err(|e| format!("cannot prepare model download staging: {e}"))?;
    let outcome = (|| {
        let mut command = client.command();
        command
            .arg("download")
            .arg(&request.repository)
            .args(&request.files)
            .args([
                "--repo-type",
                "model",
                "--revision",
                &request.revision,
                "--local-dir",
            ])
            .arg(staging.path())
            .stdin(Stdio::inherit())
            .stderr(Stdio::inherit());
        #[cfg(unix)]
        {
            use std::os::fd::AsFd;
            command.stdout(Stdio::from(
                std::io::stderr()
                    .as_fd()
                    .try_clone_to_owned()
                    .map_err(|e| format!("cannot route hf diagnostics: {e}"))?,
            ));
        }
        let status = crate::execution::run(&mut command)?;
        if !status.success() {
            return Err(format!("hf model download failed with {status}"));
        }
        let resource = write_bundle(
            staging.path(),
            &request.repository,
            &request.revision,
            &request.files,
            &version,
        )?;
        publish(staging.path(), &destination)?;
        Ok(ModelResourceResult {
            status: "downloaded",
            local_dir: destination.clone(),
            resource,
        })
    })();
    match outcome {
        Ok(result) => Ok(result),
        Err(error) => {
            let retained = staging.keep();
            // Ordinary diagnostics for incomplete, newly owned download files.
            // Never overwrite downloaded content or a successful bundle record.
            if !retained.join("BIOV_MODEL_RESOURCE.json").exists() {
                let note = serde_json::json!({
                    "status": "not_published", "provider": "huggingface", "repo_type": "model",
                    "repository": request.repository, "requested_revision": request.revision,
                    "selected_files": request.files, "hf_version": version,
                    "note": "Incomplete download; no successful resource verification is claimed. Use official hf to inspect/retry native files and its .cache/huggingface metadata."
                });
                if let Ok(mut file) = fs::OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(retained.join("BIOV_INCOMPLETE_MODEL_DOWNLOAD.json"))
                {
                    let _ = serde_json::to_writer_pretty(&mut file, &note);
                }
            }
            Err(format!(
                "{error}; unpublished files retained at {}. Existing destination is unchanged",
                retained.display()
            ))
        }
    }
}

// Apply the same direct-directory rule before acquisition and on offline reuse.
// Resolving a symlink ancestor only for a new download would publish a bundle
// that the identical request could not subsequently inspect or reuse.
fn validate_destination_components(path: &Path) -> Result<(), String> {
    let mut component_path = PathBuf::new();
    for component in path.components() {
        component_path.push(component.as_os_str());
        match fs::symlink_metadata(&component_path) {
            Ok(metadata) if metadata.is_dir() && !metadata.file_type().is_symlink() => (),
            Ok(_) => {
                return Err(format!(
                    "model-resource destination must contain only real directories: {}",
                    component_path.display()
                ))
            }
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => (),
            Err(error) => {
                return Err(format!(
                    "cannot inspect model-resource destination: {error}"
                ))
            }
        }
    }
    Ok(())
}

#[cfg(target_os = "linux")]
fn publish(source: &Path, destination: &Path) -> Result<(), String> {
    use std::{ffi::CString, os::unix::ffi::OsStrExt};
    let source = CString::new(source.as_os_str().as_bytes()).map_err(|_| "invalid staging path")?;
    let destination =
        CString::new(destination.as_os_str().as_bytes()).map_err(|_| "invalid destination path")?;
    // SAFETY: both C strings are valid for this call. Linux atomic no-replace
    // publication preserves a destination created by another downloader/user.
    let result = unsafe {
        libc::syscall(
            libc::SYS_renameat2,
            libc::AT_FDCWD,
            source.as_ptr(),
            libc::AT_FDCWD,
            destination.as_ptr(),
            libc::RENAME_NOREPLACE,
        )
    };
    if result != 0 {
        return Err(format!(
            "cannot publish model resource without replacement: {}",
            std::io::Error::last_os_error()
        ));
    }
    Ok(())
}

#[cfg(not(target_os = "linux"))]
fn publish(_source: &Path, _destination: &Path) -> Result<(), String> {
    Err("atomic model-resource publication is unsupported on this platform".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn official_version_output_forms_and_bounds() {
        for text in [
            "huggingface_hub version: 0.34.0",
            "version=1.33.0",
            "✓ hf version\n  version: 2.1.1\n",
            "2.1.1",
        ] {
            assert!(parse_version(text).is_ok(), "{text}");
        }
        for text in [
            "version: 0.33.9",
            "version: 3.0.0",
            "version: 2.1.1.dev0",
            "2.1.1\n1.0.0",
            "unrelated program",
        ] {
            assert!(parse_version(text).is_err(), "{text}");
        }
    }
}
