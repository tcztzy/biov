//! Portable companions for an explicitly selected slice of a model repository.
//!
//! This module does not resolve revisions, access the network, decode weights or
//! execute model code. The caller delegates acquisition to the official `hf`
//! client in a fresh staging directory, then records the successful invocation.
//! Checks establish local consistency, not authenticity or scientific quality.
//! Filesystem checks assume a trusted root without hostile concurrent writers.

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeSet,
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Component, Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

pub const MODEL_RESOURCE_RECORD: &str = "BIOV_MODEL_RESOURCE.json";
pub const MODEL_RESOURCE_README: &str = "BIOV_MODEL_RESOURCE_README.md";
pub const INCOMPLETE_MODEL_DOWNLOAD_RECORD: &str = "BIOV_INCOMPLETE_MODEL_DOWNLOAD.json";
pub const MODEL_RESOURCE_FORMAT_VERSION: u32 = 1;
pub const MAX_MODEL_RESOURCE_RECORD_BYTES: u64 = 1024 * 1024;
pub const MAX_MODEL_RESOURCE_FILES: usize = 1024;
pub const MAX_MODEL_RESOURCE_PATH_BYTES: usize = 512;
pub const MAX_MODEL_RESOURCE_SELECTION_BYTES: usize = 256 * 1024;
const MAX_README_BYTES: u64 = 32 * 1024;
const MAX_TREE_ENTRIES: usize = 65536;
const MAX_PATH_COMPONENTS: usize = 32;
const MAX_CLIENT_VERSION_BYTES: usize = 256;

/// Names and byte identities refer only to the selected files, relative to the
/// directory containing this record. Equality does not establish authenticity.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ModelResourceRecord {
    pub format_version: u32,
    pub provider: String,
    pub repository_type: String,
    pub repository: String,
    pub revision: String,
    pub revision_provenance: String,
    pub scope: String,
    pub selected_files: Vec<String>,
    pub inventory: Vec<ModelResourceFile>,
    pub acquisition: ModelResourceAcquisition,
    pub metadata: ModelResourceMetadata,
    pub verification: ModelResourceVerification,
    pub companions: ModelResourceCompanions,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ModelResourceFile {
    pub path: String,
    pub bytes: u64,
    pub sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ModelResourceAcquisition {
    pub client: String,
    pub client_version: String,
    pub method: String,
    pub client_invocation_status: String,
    pub companion_created_unix_seconds: u64,
    pub upstream_download_time: Option<String>,
    pub payload_transformation: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ModelResourceMetadata {
    pub model_meaning: String,
    pub native_metadata: String,
    pub model_code_execution: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ModelResourceVerification {
    pub algorithm: String,
    pub coverage: String,
    pub authenticity: String,
    pub scientific_quality_control: String,
    pub revision_resolution: String,
    pub excluded_manager_metadata: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ModelResourceCompanions {
    pub record: String,
    pub readme: String,
    pub readme_bytes: u64,
    pub readme_sha256: String,
}

/// Validate exact file selections before invoking `hf`. No globs, implicit
/// whole-repository downloads, branch names or abbreviated revisions are accepted.
pub fn validate_selection(
    repository: &str,
    revision: &str,
    files: &[String],
) -> Result<(), String> {
    let parts: Vec<_> = repository.split('/').collect();
    if repository.len() > 192
        || !(1..=2).contains(&parts.len())
        || parts.iter().any(|part| {
            part.is_empty()
                || part.len() > 96
                || part.starts_with(['.', '-'])
                || part.ends_with(['.', '-'])
                || part.contains("--")
                || part.contains("..")
                || part.ends_with(".git")
                || !part
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b))
        })
    {
        return Err(
            "model repository must be a bounded Hugging Face name or namespace/name".into(),
        );
    }
    if revision.len() != 40 || !revision.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err(
            "model revision must be an exact caller-supplied full 40-hex Git commit".into(),
        );
    }
    if files.is_empty() || files.len() > MAX_MODEL_RESOURCE_FILES {
        return Err(format!(
            "model selection must contain 1..={MAX_MODEL_RESOURCE_FILES} explicit files"
        ));
    }
    // A conservative aggregate path budget lets callers reject an oversized
    // record before any expensive download. Paths appear twice in the record.
    if files.iter().map(String::len).sum::<usize>() > MAX_MODEL_RESOURCE_SELECTION_BYTES {
        return Err("model selection exceeds the 256 KiB aggregate path byte bound".into());
    }
    let mut names = BTreeSet::new();
    for name in files {
        validate_relative_file(name)?;
        if !names.insert(name.as_str()) {
            return Err(format!("duplicate model file selection: {name}"));
        }
    }
    for name in &names {
        let components: Vec<_> = name.split('/').collect();
        for end in 1..components.len() {
            if names.contains(components[..end].join("/").as_str()) {
                return Err(format!(
                    "model file selection has a file/directory collision: {name}"
                ));
            }
        }
    }
    Ok(())
}

fn validate_relative_file(name: &str) -> Result<(), String> {
    let components: Vec<_> = name.split('/').collect();
    if name.is_empty()
        || name.len() > MAX_MODEL_RESOURCE_PATH_BYTES
        || components.len() > MAX_PATH_COMPONENTS
        || components
            .iter()
            .any(|part| part.is_empty() || *part == "." || *part == "..")
        || name
            .bytes()
            .any(|b| b.is_ascii_control() || b"\\:*?[]{}<>\"|".contains(&b))
        || name.starts_with('-')
        || components[0].eq_ignore_ascii_case(".cache")
        || components.iter().any(|part| {
            part.eq_ignore_ascii_case(MODEL_RESOURCE_RECORD)
                || part.eq_ignore_ascii_case(MODEL_RESOURCE_README)
                || part.eq_ignore_ascii_case(INCOMPLETE_MODEL_DOWNLOAD_RECORD)
        })
    {
        return Err(format!(
            "unsafe, reserved or non-explicit model file path: {name:?}"
        ));
    }
    if Path::new(name)
        .components()
        .any(|c| !matches!(c, Component::Normal(_)))
    {
        return Err(format!("model file path must be relative: {name:?}"));
    }
    Ok(())
}

fn validate_client_version(version: &str) -> Result<(), String> {
    if version.is_empty()
        || version.len() > MAX_CLIENT_VERSION_BYTES
        || version.trim() != version
        || version.chars().any(char::is_control)
    {
        return Err("hf client version must be a nonempty bounded version string".into());
    }
    Ok(())
}

/// Add companions without changing payloads or overwriting existing companions.
///
/// Precondition: the caller observed the official `hf download` command succeed
/// in a fresh staging root with this model repository, exact revision and explicit
/// file list. `hf_version` is the actual invoked client's reported version.
/// These acquisition facts are the calling manager's evidence, not an independent
/// upstream audit. A failure may leave a newly created README in the staging root;
/// the record is published last, and callers must not publish failed staging roots.
pub fn write_bundle(
    root: &Path,
    repository: &str,
    revision: &str,
    files: &[String],
    hf_version: &str,
) -> Result<ModelResourceRecord, String> {
    validate_selection(repository, revision, files)?;
    validate_client_version(hf_version)?;
    let root = checked_root(root)?;
    for name in [MODEL_RESOURCE_RECORD, MODEL_RESOURCE_README] {
        match fs::symlink_metadata(root.join(name)) {
            Ok(_) => return Err(format!("refusing to overwrite model companion: {name}")),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(format!("cannot check model companion {name}: {e}")),
        }
    }
    let mut selected_files = files.to_vec();
    selected_files.sort();
    check_payload_tree(&root, &selected_files, false)?;
    let inventory = selected_files
        .iter()
        .map(|name| hash_payload(&root, name))
        .collect::<Result<Vec<_>, _>>()?;
    let created = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| format!("cannot record companion creation time: {e}"))?
        .as_secs();
    let record = make_record(
        repository,
        revision,
        selected_files,
        inventory,
        hf_version,
        created,
    );
    let bytes = serde_json::to_vec_pretty(&record).map_err(|e| e.to_string())?;
    if bytes.len() as u64 > MAX_MODEL_RESOURCE_RECORD_BYTES {
        return Err("model resource record exceeds the 1 MiB bound".into());
    }
    let readme = prepared_companion(&root, README.as_bytes())?;
    let manifest = prepared_companion(&root, &bytes)?;
    readme
        .persist_noclobber(root.join(MODEL_RESOURCE_README))
        .map_err(|e| format!("cannot publish model README without replacing a file: {e}"))?;
    manifest
        .persist_noclobber(root.join(MODEL_RESOURCE_RECORD))
        .map_err(|e| format!("cannot publish model record without replacing a file: {e}"))?;
    Ok(record)
}

/// Read and verify complete selected-file bytes offline. No original directory,
/// cache database, Hugging Face installation or network access is required.
pub fn inspect_bundle(root: &Path) -> Result<ModelResourceRecord, String> {
    let root = checked_root(root)?;
    let bytes = read_bounded_regular(
        &root,
        MODEL_RESOURCE_RECORD,
        MAX_MODEL_RESOURCE_RECORD_BYTES,
    )?;
    let record: ModelResourceRecord = serde_json::from_slice(&bytes)
        .map_err(|e| format!("invalid model resource record: {e}"))?;
    validate_record(&record)?;
    let readme = read_bounded_regular(&root, MODEL_RESOURCE_README, MAX_README_BYTES)?;
    if readme != README.as_bytes() {
        return Err("model resource README does not match the supported portable format".into());
    }
    check_payload_tree(&root, &record.selected_files, true)?;
    for expected in &record.inventory {
        let actual = hash_payload(&root, &expected.path)?;
        if actual != *expected {
            return Err(format!(
                "model payload size or SHA-256 mismatch: {}",
                expected.path
            ));
        }
    }
    Ok(record)
}

fn make_record(
    repository: &str,
    revision: &str,
    selected_files: Vec<String>,
    inventory: Vec<ModelResourceFile>,
    hf_version: &str,
    created: u64,
) -> ModelResourceRecord {
    ModelResourceRecord {
        format_version: MODEL_RESOURCE_FORMAT_VERSION,
        provider: "hugging_face".into(),
        repository_type: "model".into(),
        repository: repository.into(),
        revision: revision.into(),
        revision_provenance: "caller_supplied_full_git_commit_passed_to_hf".into(),
        scope: "selected_files".into(),
        selected_files,
        inventory,
        acquisition: ModelResourceAcquisition {
            client: "hf".into(),
            client_version: hf_version.into(),
            method: "hf_download_explicit_files_exact_revision_fresh_local_dir".into(),
            client_invocation_status: "success_observed_by_calling_manager".into(),
            companion_created_unix_seconds: created,
            upstream_download_time: None,
            payload_transformation: "none_by_biov".into(),
        },
        metadata: ModelResourceMetadata {
            model_meaning: "not_interpreted".into(),
            native_metadata: "authoritative_when_present_in_selected_files".into(),
            model_code_execution: "not_performed".into(),
        },
        verification: ModelResourceVerification {
            algorithm: "sha256".into(),
            coverage: "complete_bytes_of_each_selected_file".into(),
            authenticity: "not_established".into(),
            scientific_quality_control: "not_performed".into(),
            revision_resolution: "not_independently_performed".into(),
            excluded_manager_metadata: vec![".cache/huggingface".into()],
        },
        companions: ModelResourceCompanions {
            record: MODEL_RESOURCE_RECORD.into(),
            readme: MODEL_RESOURCE_README.into(),
            readme_bytes: README.len() as u64,
            readme_sha256: format!("{:x}", Sha256::digest(README.as_bytes())),
        },
    }
}

fn validate_record(record: &ModelResourceRecord) -> Result<(), String> {
    validate_selection(&record.repository, &record.revision, &record.selected_files)?;
    validate_client_version(&record.acquisition.client_version)?;
    if record.inventory.len() != record.selected_files.len() {
        return Err("model inventory must cover exactly the selected files".into());
    }
    for (index, file) in record.inventory.iter().enumerate() {
        if file.path != record.selected_files[index]
            || file.sha256.len() != 64
            || !file
                .sha256
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            || (index > 0 && record.selected_files[index - 1] >= record.selected_files[index])
        {
            return Err(
                "model inventory must be sorted, exact and contain lowercase SHA-256 values".into(),
            );
        }
    }
    let expected = make_record(
        &record.repository,
        &record.revision,
        record.selected_files.clone(),
        record.inventory.clone(),
        &record.acquisition.client_version,
        record.acquisition.companion_created_unix_seconds,
    );
    if *record != expected {
        return Err(
            "unsupported or inconsistent model record schema, scope or acquisition facts".into(),
        );
    }
    Ok(())
}

fn checked_root(root: &Path) -> Result<PathBuf, String> {
    let absolute = if root.is_absolute() {
        root.to_path_buf()
    } else {
        std::env::current_dir()
            .map_err(|e| e.to_string())?
            .join(root)
    };
    let mut part = PathBuf::new();
    for component in absolute.components() {
        part.push(component.as_os_str());
        let metadata = fs::symlink_metadata(&part)
            .map_err(|e| format!("cannot inspect model root {}: {e}", part.display()))?;
        if metadata.file_type().is_symlink() || !metadata.is_dir() {
            return Err(format!(
                "model root must contain only real directories: {}",
                part.display()
            ));
        }
    }
    fs::canonicalize(&absolute).map_err(|e| format!("cannot resolve model root: {e}"))
}

fn checked_file(root: &Path, name: &str) -> Result<PathBuf, String> {
    let mut path = root.to_path_buf();
    let parts: Vec<_> = name.split('/').collect();
    for (index, part) in parts.iter().enumerate() {
        path.push(part);
        let metadata = fs::symlink_metadata(&path)
            .map_err(|e| format!("cannot inspect model file {name}: {e}"))?;
        if metadata.file_type().is_symlink()
            || (index + 1 == parts.len() && !metadata.is_file())
            || (index + 1 < parts.len() && !metadata.is_dir())
        {
            return Err(format!(
                "model payload must be a regular file without symlink components: {name}"
            ));
        }
    }
    let resolved =
        fs::canonicalize(&path).map_err(|e| format!("cannot resolve model file {name}: {e}"))?;
    if !resolved.starts_with(root) {
        return Err(format!("model file escapes its root: {name}"));
    }
    Ok(path)
}

fn open_regular(path: &Path) -> Result<File, String> {
    let mut options = OpenOptions::new();
    options.read(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        options.custom_flags(libc::O_NOFOLLOW | libc::O_NONBLOCK);
    }
    let file = options
        .open(path)
        .map_err(|e| format!("cannot open model file {}: {e}", path.display()))?;
    if !file.metadata().map_err(|e| e.to_string())?.is_file() {
        return Err(format!("model file is not regular: {}", path.display()));
    }
    Ok(file)
}

fn hash_payload(root: &Path, name: &str) -> Result<ModelResourceFile, String> {
    let path = checked_file(root, name)?;
    let mut file = open_regular(&path)?;
    let before = file.metadata().map_err(|e| e.to_string())?.len();
    let mut sha = Sha256::new();
    let mut bytes = 0_u64;
    let mut buffer = [0_u8; 128 * 1024];
    loop {
        let n = file
            .read(&mut buffer)
            .map_err(|e| format!("cannot read model payload {name}: {e}"))?;
        if n == 0 {
            break;
        }
        bytes = bytes
            .checked_add(n as u64)
            .ok_or("model payload byte count overflow")?;
        sha.update(&buffer[..n]);
    }
    if before != bytes || file.metadata().map_err(|e| e.to_string())?.len() != bytes {
        return Err(format!(
            "model payload changed length during verification: {name}"
        ));
    }
    Ok(ModelResourceFile {
        path: name.into(),
        bytes,
        sha256: format!("{:x}", sha.finalize()),
    })
}

fn read_bounded_regular(root: &Path, name: &str, limit: u64) -> Result<Vec<u8>, String> {
    let path = checked_file(root, name)?;
    let mut file = open_regular(&path)?;
    if file.metadata().map_err(|e| e.to_string())?.len() > limit {
        return Err(format!(
            "model companion exceeds bounded byte limit: {name}"
        ));
    }
    let mut bytes = Vec::new();
    Read::by_ref(&mut file)
        .take(limit + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    if bytes.len() as u64 > limit {
        return Err(format!(
            "model companion exceeds bounded byte limit: {name}"
        ));
    }
    Ok(bytes)
}

fn check_payload_tree(
    root: &Path,
    selected: &[String],
    allow_companions: bool,
) -> Result<(), String> {
    let expected: BTreeSet<_> = selected.iter().map(String::as_str).collect();
    let mut seen = BTreeSet::new();
    let mut pending = vec![(root.to_path_buf(), 0_usize)];
    let mut entries = 0_usize;
    while let Some((directory, depth)) = pending.pop() {
        for entry in
            fs::read_dir(&directory).map_err(|e| format!("cannot list model bundle: {e}"))?
        {
            let entry = entry.map_err(|e| e.to_string())?;
            entries += 1;
            if entries > MAX_TREE_ENTRIES {
                return Err("model bundle directory traversal exceeds its entry bound".into());
            }
            let path = entry.path();
            let relative = path.strip_prefix(root).map_err(|e| e.to_string())?;
            let name = relative
                .to_str()
                .ok_or("model bundle paths must be UTF-8")?
                .replace(std::path::MAIN_SEPARATOR, "/");
            let metadata = fs::symlink_metadata(&path).map_err(|e| e.to_string())?;
            if metadata.file_type().is_symlink() {
                return Err(format!("model bundle contains a symlink: {name}"));
            }
            if name == ".cache/huggingface" {
                if !metadata.is_dir() {
                    return Err("Hugging Face manager metadata must be a directory".into());
                }
                // The client's private metadata is excluded, never followed or
                // required for reading, and is not part of the payload inventory.
                continue;
            }
            if name.starts_with(".cache/") {
                return Err(format!(
                    "unexpected file outside the excluded hf metadata subtree: {name}"
                ));
            }
            if metadata.is_dir() {
                if depth >= MAX_PATH_COMPONENTS {
                    return Err("model bundle directory depth exceeds its bound".into());
                }
                pending.push((path, depth + 1));
            } else if metadata.is_file() {
                if allow_companions
                    && [MODEL_RESOURCE_RECORD, MODEL_RESOURCE_README].contains(&name.as_str())
                {
                    continue;
                }
                if !expected.contains(name.as_str()) {
                    return Err(format!(
                        "model bundle contains an unrecorded payload: {name}"
                    ));
                }
                seen.insert(name);
            } else {
                return Err(format!(
                    "model bundle contains a nonregular payload: {name}"
                ));
            }
        }
    }
    if seen.len() != expected.len() || seen.iter().any(|name| !expected.contains(name.as_str())) {
        return Err("model bundle is missing one or more selected regular payload files".into());
    }
    Ok(())
}

fn prepared_companion(root: &Path, bytes: &[u8]) -> Result<tempfile::NamedTempFile, String> {
    let mut file = tempfile::Builder::new()
        .prefix(".biov-model-companion-")
        .tempfile_in(root)
        .map_err(|e| e.to_string())?;
    file.write_all(bytes).map_err(|e| e.to_string())?;
    file.as_file().sync_all().map_err(|e| e.to_string())?;
    Ok(file)
}

// Format-1 companion bytes are immutable. Edit user-facing guides separately;
// changing this template requires deliberate versioning and old-format support.
const README: &str = r#"# Portable selected model files

This directory retains the selected provider-native files and relative layout.
BIOV_MODEL_RESOURCE.json is the machine-readable entry point (format_version 1).
BIOV_MODEL_RESOURCE_README.md is this guide. Keep both companions and every
inventory file together when copying or moving this directory. BioV, the original
directory, Hugging Face, a cache catalog and network access are not required to
read or verify it. The optional .cache/huggingface subtree is private client
metadata: it is excluded from verification and is not required for reuse.

## What is known

The record identifies a Hugging Face model repository, the exact full 40-hex Git
revision supplied by the caller, and the explicit selected files. scope is always
selected_files. This is not a claim that the complete model, all weight shards,
configuration, tokenizer, license, documentation or execution dependencies are
present. Consult selected provider-native documentation and metadata when present.
No model architecture, training data, reference assembly, biological meaning,
license or scientific fitness is inferred from file names or repository names.

inventory lists every selected relative path, complete byte count and SHA-256.
No payload is converted or decoded by BioV, and no downloaded model code is run.
Use the appropriate ordinary reader for a native format only after establishing
what that format contains. In particular, this example never unpickles weights
or imports Python from the bundle.

The acquisition facts report the actual hf client/version and a successful
explicit-file, exact-revision download invocation observed by the calling manager
in a fresh staging directory. The companion creation time is a local recording
time, not an upstream publication or original download time; upstream_download_time
is unknown (null). The Git revision is caller supplied and was passed to hf; BioV
does not independently resolve or authenticate the provider's returned revision.
The client version and format_version are software/schema versions, not biological
versions. Checksums establish consistency with this supplied record, not producer
authenticity, verified source provenance or scientific quality control. A changed
record can describe changed bytes; signatures and source authentication are outside
this contract. Filesystem checks assume no hostile concurrent filesystem writer.

## Independent complete-byte verification and summary

The following uses only the Python standard library. It validates the bounded
record schema, the README byte identity, exact selected-file inventory, safe relative paths and complete
payload bytes, then summarizes all inventory records by native filename suffix.
It never contacts the network and never loads a model. It rejects unrecorded
payloads, symlink paths, changed/missing files and malformed companions. The private
hf metadata subtree is ignored. The reported suffix groups have no scientific
meaning. Run from any directory, including after removing BioV:

```sh
python3 -I -S - /absolute/path/to/moved-bundle <<'PY'
import hashlib
import json
import os
from pathlib import Path
import sys

RECORD = 'BIOV_MODEL_RESOURCE.json'
README = 'BIOV_MODEL_RESOURCE_README.md'
MAX_RECORD = 1024 * 1024
MAX_FILES = 1024
MAX_PATH = 512
MAX_SELECTION_BYTES = 256 * 1024
MAX_DEPTH = 32
MAX_ENTRIES = 65536

def require(condition, message):
    if not condition:
        raise ValueError(message)

root = Path(os.path.abspath(sys.argv[1]))
for component in (root, *root.parents):
    require(not component.is_symlink() and component.is_dir(), 'unsafe root')
root = root.resolve()

def regular(relative):
    path = root
    parts = relative.split('/')
    for index, part in enumerate(parts):
        path = path / part
        require(not path.is_symlink(), 'symlink: ' + relative)
        require(path.is_file() if index == len(parts) - 1 else path.is_dir(),
                'nonregular or missing path: ' + relative)
    require(root in path.resolve().parents, 'path escape')
    return path

manifest_path = regular(RECORD)
with manifest_path.open('rb') as stream:
    raw = stream.read(MAX_RECORD + 1)
require(len(raw) <= MAX_RECORD, 'record byte bound')
# Reject duplicate JSON keys instead of silently taking the last value.
def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'duplicate JSON field')
        result[key] = value
    return result
record = json.loads(raw.decode('utf-8'), object_pairs_hook=unique_object)
require(type(record) is dict, 'record must be an object')
repository = record['repository']
revision = record['revision']
files = record['selected_files']
require(type(repository) is str and len(repository.encode('utf-8')) <= 192,
        'repository bound')
parts = repository.split('/')
require(1 <= len(parts) <= 2, 'repository syntax')
for part in parts:
    require(1 <= len(part) <= 96 and part[0] not in '.-' and part[-1] not in '.-'
            and '--' not in part and '..' not in part and not part.endswith('.git')
            and all(c.isascii() and (c.isalnum() or c in '._-') for c in part),
            'repository syntax')
require(type(revision) is str and len(revision) == 40
        and all(c in '0123456789abcdefABCDEF' for c in revision), 'full commit required')
require(type(files) is list and 1 <= len(files) <= MAX_FILES, 'selection bound')
for name in files:
    require(type(name) is str and 1 <= len(name.encode('utf-8')) <= MAX_PATH,
            'file path bound')
    parts = name.split('/')
    require(len(parts) <= MAX_DEPTH and all(p not in ('', '.', '..') for p in parts)
            and not name.startswith('-')
            and parts[0].lower() != '.cache'
            and all(p.lower() not in (RECORD.lower(), README.lower(),
                                     'biov_incomplete_model_download.json') for p in parts)
            and all(ord(c) >= 32 and ord(c) != 127 and c not in '\\:*?[]{}<>"|'
                    for c in name), 'unsafe or reserved file selection')
require(sum(len(name.encode('utf-8')) for name in files) <= MAX_SELECTION_BYTES,
        'aggregate selected path byte bound')
require(files == sorted(set(files)), 'selection must be sorted and unique')
for name in files:
    parts = name.split('/')
    require(not any('/'.join(parts[:i]) in files for i in range(1, len(parts))),
            'file/directory collision')
inventory = record['inventory']
require(type(inventory) is list and len(inventory) == len(files), 'inventory coverage')
for name, item in zip(files, inventory):
    require(type(item) is dict and set(item) == {'path', 'bytes', 'sha256'}
            and item['path'] == name and type(item['bytes']) is int
            and 0 <= item['bytes'] <= 2**64 - 1
            and type(item['sha256']) is str and len(item['sha256']) == 64
            and all(c in '0123456789abcdef' for c in item['sha256']), 'inventory schema')
acquisition = record['acquisition']
version = acquisition['client_version']
created = acquisition['companion_created_unix_seconds']
require(type(version) is str and 1 <= len(version.encode('utf-8')) <= 256
        and version == version.strip() and all(not (ord(c) < 32 or 127 <= ord(c) <= 159)
                                                for c in version), 'client version')
require(type(created) is int and 0 <= created <= 2**64 - 1, 'creation time')
companions = record['companions']
readme_bytes = companions['readme_bytes']
readme_sha256 = companions['readme_sha256']
require(type(readme_bytes) is int and 1 <= readme_bytes <= 32 * 1024
        and type(readme_sha256) is str and len(readme_sha256) == 64
        and all(c in '0123456789abcdef' for c in readme_sha256), 'README identity')
expected = {
    'format_version': 1, 'provider': 'hugging_face', 'repository_type': 'model',
    'repository': repository, 'revision': revision,
    'revision_provenance': 'caller_supplied_full_git_commit_passed_to_hf',
    'scope': 'selected_files', 'selected_files': files, 'inventory': inventory,
    'acquisition': {
        'client': 'hf', 'client_version': version,
        'method': 'hf_download_explicit_files_exact_revision_fresh_local_dir',
        'client_invocation_status': 'success_observed_by_calling_manager',
        'companion_created_unix_seconds': created, 'upstream_download_time': None,
        'payload_transformation': 'none_by_biov',
    },
    'metadata': {
        'model_meaning': 'not_interpreted',
        'native_metadata': 'authoritative_when_present_in_selected_files',
        'model_code_execution': 'not_performed',
    },
    'verification': {
        'algorithm': 'sha256', 'coverage': 'complete_bytes_of_each_selected_file',
        'authenticity': 'not_established', 'scientific_quality_control': 'not_performed',
        'revision_resolution': 'not_independently_performed',
        'excluded_manager_metadata': ['.cache/huggingface'],
    },
    'companions': {'record': RECORD, 'readme': README,
                   'readme_bytes': readme_bytes, 'readme_sha256': readme_sha256},
}
require(type(record['format_version']) is int and record == expected, 'record schema or scope')
with regular(README).open('rb') as stream:
    readme_raw = stream.read(32 * 1024 + 1)
require(len(readme_raw) == readme_bytes
        and hashlib.sha256(readme_raw).hexdigest() == readme_sha256, 'README identity mismatch')

seen = set()
pending = [(root, 0)]
entries = 0
while pending:
    directory, depth = pending.pop()
    for path in directory.iterdir():
        entries += 1
        require(entries <= MAX_ENTRIES, 'directory entry bound')
        name = path.relative_to(root).as_posix()
        require(not path.is_symlink(), 'symlink: ' + name)
        if name == '.cache/huggingface':
            require(path.is_dir(), 'hf metadata must be a directory')
            continue
        require(not name.startswith('.cache/'), 'unexpected private metadata')
        if path.is_dir():
            require(depth < MAX_DEPTH, 'directory depth bound')
            pending.append((path, depth + 1))
        else:
            require(path.is_file(), 'nonregular file: ' + name)
            if name not in (RECORD, README):
                require(name in files, 'unrecorded payload: ' + name)
                seen.add(name)
require(seen == set(files), 'missing selected files')

groups = {}
total_bytes = 0
for item in inventory:
    digest = hashlib.sha256()
    count = 0
    with regular(item['path']).open('rb') as stream:
        while True:
            block = stream.read(128 * 1024)
            if not block:
                break
            count += len(block)
            digest.update(block)
    require(count == item['bytes'] and digest.hexdigest() == item['sha256'],
            'byte count or SHA-256 mismatch: ' + item['path'])
    total_bytes += count
    suffix = Path(item['path']).suffix or '(no suffix)'
    group = groups.setdefault(suffix, {'files': 0, 'bytes': 0})
    group['files'] += 1
    group['bytes'] += count
print(json.dumps({'scope': record['scope'], 'verified_files': len(inventory),
                  'complete_payload_bytes': total_bytes,
                  'native_filename_suffixes': groups}, sort_keys=True))
PY
```

Copy the complete shell block, or save the Python between the command and PY
as verify_model_resource.py and run python3 -I -S verify_model_resource.py BUNDLE.
The README and JSON record are ordinary text; no BioV-specific reader is needed.
Bounds: 1 MiB JSON record, 1,024 selected files, 512 UTF-8 bytes per relative path,
256 KiB of aggregate selected path bytes,
32 path components, 65,536 traversed directory entries. Hashing streams whole files
with a 128 KiB buffer; there is no whole-weight allocation or payload size cap.
"#;

#[cfg(test)]
mod tests {
    use super::*;
    use std::process::Command;

    const REVISION: &str = "0123456789abcdef0123456789abcdef01234567";

    fn fixture() -> (tempfile::TempDir, Vec<String>) {
        let root = tempfile::tempdir().unwrap();
        fs::create_dir(root.path().join("weights")).unwrap();
        fs::write(root.path().join("config.json"), b"{\"example\":true}\n").unwrap();
        // More than one hashing buffer, including zeros and nontext bytes.
        fs::write(
            root.path().join("weights/model.safetensors"),
            vec![0xa5; 300_001],
        )
        .unwrap();
        let files = vec!["weights/model.safetensors".into(), "config.json".into()];
        (root, files)
    }

    fn write_fixture(root: &Path, files: &[String]) -> ModelResourceRecord {
        write_bundle(root, "example/selected-model", REVISION, files, "1.2.3").unwrap()
    }

    fn change_record(root: &Path, change: impl FnOnce(&mut serde_json::Value)) {
        let path = root.join(MODEL_RESOURCE_RECORD);
        let mut value: serde_json::Value =
            serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        change(&mut value);
        fs::write(path, serde_json::to_vec(&value).unwrap()).unwrap();
    }

    #[test]
    fn validates_only_exact_bounded_selections() {
        assert!(validate_selection(
            "owner/model",
            REVISION,
            &["weights/model.safetensors".into()]
        )
        .is_ok());
        assert!(
            validate_selection("model", &REVISION.to_uppercase(), &["README.md".into()]).is_ok()
        );
        for revision in [
            "main",
            "v1.0",
            "01234567",
            "g123456789abcdef0123456789abcdef01234567",
        ] {
            assert!(validate_selection("owner/model", revision, &["config.json".into()]).is_err());
        }
        for repository in [
            "",
            "/model",
            "../model",
            "owner/model/more",
            "owner/model.git",
            "a/b--c",
            "a/b c",
        ] {
            assert!(validate_selection(repository, REVISION, &["config.json".into()]).is_err());
        }
        for file in [
            "",
            "/absolute",
            "../outside",
            "a/../b",
            "./config",
            "a//b",
            "a/",
            "C:/model",
            "a\\b",
            "*.bin",
            "a[1].bin",
            ".cache/huggingface/x",
            ".cache",
            MODEL_RESOURCE_RECORD,
            MODEL_RESOURCE_README,
            "a/BIOV_MODEL_RESOURCE.json",
            "-option",
            "a\n.bin",
        ] {
            assert!(
                validate_selection("a/b", REVISION, &[file.into()]).is_err(),
                "accepted {file:?}"
            );
        }
        assert!(validate_selection("a/b", REVISION, &[]).is_err());
        assert!(validate_selection("a/b", REVISION, &["a".into(), "a".into()]).is_err());
        assert!(validate_selection("a/b", REVISION, &["a".into(), "a/b".into()]).is_err());
        assert!(validate_selection(
            "a/b",
            REVISION,
            &["a".repeat(MAX_MODEL_RESOURCE_PATH_BYTES + 1)]
        )
        .is_err());
        let many: Vec<_> = (0..=MAX_MODEL_RESOURCE_FILES)
            .map(|i| format!("file{i}"))
            .collect();
        assert!(validate_selection("a/b", REVISION, &many).is_err());
        let large_names: Vec<_> = (0..MAX_MODEL_RESOURCE_FILES)
            .map(|i| format!("{}-{i}", "a".repeat(300)))
            .collect();
        assert!(validate_selection("a/b", REVISION, &large_names)
            .unwrap_err()
            .contains("aggregate path byte bound"));
    }

    #[test]
    fn preserves_native_bytes_and_has_accurate_selected_only_facts() {
        let (root, files) = fixture();
        fs::create_dir_all(root.path().join(".cache/huggingface/download")).unwrap();
        fs::write(
            root.path()
                .join(".cache/huggingface/download/config.metadata"),
            b"client state",
        )
        .unwrap();
        let before = fs::read(root.path().join(&files[0])).unwrap();
        let record = write_fixture(root.path(), &files);
        assert_eq!(record.scope, "selected_files");
        assert_eq!(
            record.selected_files,
            vec!["config.json", "weights/model.safetensors"]
        );
        assert_eq!(record.inventory.len(), 2);
        assert_eq!(record.revision, REVISION);
        assert_eq!(record.acquisition.client_version, "1.2.3");
        assert_eq!(record.acquisition.upstream_download_time, None);
        assert_eq!(record.verification.authenticity, "not_established");
        assert_eq!(record.metadata.model_meaning, "not_interpreted");
        assert_eq!(fs::read(root.path().join(&files[0])).unwrap(), before);
        assert_eq!(inspect_bundle(root.path()).unwrap(), record);
        let mut encoded = serde_json::to_value(&record).unwrap();
        encoded["scope"] = "complete_model".into();
        let invalid: ModelResourceRecord = serde_json::from_value(encoded).unwrap();
        assert!(validate_record(&invalid).is_err());
    }

    #[test]
    fn rejects_changed_missing_and_unrecorded_payloads() {
        let (root, files) = fixture();
        write_fixture(root.path(), &files);
        let payload = root.path().join(&files[0]);
        let original = fs::read(&payload).unwrap();
        let mut changed = original.clone();
        changed[250_000] ^= 1; // Same size; whole-file hashing catches a late change.
        fs::write(&payload, changed).unwrap();
        assert!(inspect_bundle(root.path())
            .unwrap_err()
            .contains("SHA-256 mismatch"));
        fs::write(&payload, &original).unwrap();
        fs::write(
            root.path().join("extra.py"),
            b"raise RuntimeError('never execute')",
        )
        .unwrap();
        assert!(inspect_bundle(root.path())
            .unwrap_err()
            .contains("unrecorded payload"));
        fs::remove_file(root.path().join("extra.py")).unwrap();
        fs::remove_file(&payload).unwrap();
        assert!(inspect_bundle(root.path()).unwrap_err().contains("missing"));
    }

    #[test]
    fn rejects_existing_companions_without_overwriting_anything() {
        for name in [MODEL_RESOURCE_RECORD, MODEL_RESOURCE_README] {
            let (root, files) = fixture();
            fs::write(root.path().join(name), b"retain this existing companion").unwrap();
            assert!(write_bundle(root.path(), "a/b", REVISION, &files, "1.2.3").is_err());
            assert_eq!(
                fs::read(root.path().join(name)).unwrap(),
                b"retain this existing companion"
            );
            assert_eq!(
                fs::read(root.path().join(&files[0])).unwrap().len(),
                300_001
            );
        }
        let (root, files) = fixture();
        let first = write_fixture(root.path(), &files);
        assert!(write_bundle(root.path(), "a/b", REVISION, &files, "1.2.3").is_err());
        assert_eq!(inspect_bundle(root.path()).unwrap(), first);
    }

    #[test]
    fn rejects_incomplete_staging_and_invalid_client_version() {
        let (root, files) = fixture();
        fs::write(root.path().join("unselected.bin"), b"extra").unwrap();
        assert!(write_bundle(root.path(), "a/b", REVISION, &files, "1.2.3").is_err());
        fs::remove_file(root.path().join("unselected.bin")).unwrap();
        for version in ["", " 1.0", "1.0\n", "\u{85}"] {
            assert!(write_bundle(root.path(), "a/b", REVISION, &files, version).is_err());
        }
        fs::remove_file(root.path().join(&files[1])).unwrap();
        assert!(write_bundle(root.path(), "a/b", REVISION, &files, "1.2.3").is_err());
        assert!(!root.path().join(MODEL_RESOURCE_RECORD).exists());
    }

    #[test]
    fn rejects_malformed_records_scope_inventory_and_bounds() {
        for alteration in 0..7 {
            let (root, files) = fixture();
            write_fixture(root.path(), &files);
            change_record(root.path(), |value| match alteration {
                0 => value["format_version"] = 99.into(),
                1 => value["scope"] = "complete_model".into(),
                2 => value["inventory"][0]["path"] = "../escaped".into(),
                3 => value["inventory"][0]["sha256"] = "f".into(),
                4 => value["inventory"] = serde_json::json!([]),
                5 => value["acquisition"]["upstream_download_time"] = "invented".into(),
                _ => value["extra"] = true.into(),
            });
            assert!(inspect_bundle(root.path()).is_err());
        }
        let (root, files) = fixture();
        write_fixture(root.path(), &files);
        fs::write(
            root.path().join(MODEL_RESOURCE_RECORD),
            vec![b' '; MAX_MODEL_RESOURCE_RECORD_BYTES as usize + 1],
        )
        .unwrap();
        assert!(inspect_bundle(root.path())
            .unwrap_err()
            .contains("bounded byte limit"));
        fs::write(root.path().join(MODEL_RESOURCE_RECORD), b"{broken JSON").unwrap();
        assert!(inspect_bundle(root.path())
            .unwrap_err()
            .contains("invalid model resource record"));
    }

    #[test]
    fn rejects_modified_readme_and_duplicate_json_fields() {
        let (root, files) = fixture();
        write_fixture(root.path(), &files);
        let readme_path = root.path().join(MODEL_RESOURCE_README);
        let original_readme = fs::read(&readme_path).unwrap();
        fs::write(&readme_path, b"changed companion").unwrap();
        assert!(inspect_bundle(root.path()).unwrap_err().contains("README"));
        fs::write(readme_path, original_readme).unwrap();
        let record_path = root.path().join(MODEL_RESOURCE_RECORD);
        let json = fs::read_to_string(&record_path).unwrap();
        let duplicate = json.replacen('{', "{\"format_version\":1,", 1);
        fs::write(record_path, duplicate).unwrap();
        assert!(inspect_bundle(root.path())
            .unwrap_err()
            .contains("duplicate field"));
    }

    #[cfg(unix)]
    #[test]
    fn rejects_payload_and_directory_symlinks_and_nonregular_files() {
        use std::os::unix::fs::symlink;
        let outside = tempfile::tempdir().unwrap();
        fs::write(outside.path().join("weights.bin"), b"outside").unwrap();
        let root = tempfile::tempdir().unwrap();
        symlink(
            outside.path().join("weights.bin"),
            root.path().join("weights.bin"),
        )
        .unwrap();
        assert!(write_bundle(
            root.path(),
            "a/b",
            REVISION,
            &["weights.bin".into()],
            "1.2.3"
        )
        .is_err());
        fs::remove_file(root.path().join("weights.bin")).unwrap();
        symlink(outside.path(), root.path().join("link")).unwrap();
        assert!(write_bundle(
            root.path(),
            "a/b",
            REVISION,
            &["link/weights.bin".into()],
            "1.2.3"
        )
        .is_err());
        let alias = outside.path().join("root-link");
        symlink(root.path(), &alias).unwrap();
        assert!(write_bundle(
            &alias,
            "a/b",
            REVISION,
            &["link/weights.bin".into()],
            "1.2.3"
        )
        .is_err());
        fs::remove_file(root.path().join("link")).unwrap();
        let fifo = root.path().join("pipe");
        use std::os::unix::ffi::OsStrExt;
        let fifo_name = std::ffi::CString::new(fifo.as_os_str().as_bytes()).unwrap();
        assert_eq!(unsafe { libc::mkfifo(fifo_name.as_ptr(), 0o600) }, 0);
        assert!(write_bundle(root.path(), "a/b", REVISION, &["pipe".into()], "1.2.3").is_err());
    }

    #[test]
    fn moved_bundle_has_independent_offline_complete_byte_reader() {
        let python = match Command::new("python3")
            .args(["-I", "-S", "-c", "import sys; print(sys.executable)"])
            .output()
        {
            Ok(output) if output.status.success() => {
                String::from_utf8(output.stdout).unwrap().trim().to_owned()
            }
            _ => {
                eprintln!("python3 unavailable; independent-reader acceptance test skipped");
                return;
            }
        };
        let (original, files) = fixture();
        let record = write_fixture(original.path(), &files);
        let moved_parent = tempfile::tempdir().unwrap();
        let moved = moved_parent.path().join("moved-native-bundle");
        fs::rename(original.path(), &moved).unwrap();
        drop(original);
        assert_eq!(inspect_bundle(&moved).unwrap(), record);
        let readme = fs::read_to_string(moved.join(MODEL_RESOURCE_README)).unwrap();
        let script = readme
            .split("<<'PY'\n")
            .nth(1)
            .unwrap()
            .split("\nPY\n")
            .next()
            .unwrap();
        let isolated_script = format!(
            "{script}\nimport importlib.util\nassert importlib.util.find_spec('biov') is None\n"
        );
        let run = || {
            Command::new(&python)
                .args(["-I", "-S", "-c", &isolated_script])
                .arg(&moved)
                .env_clear()
                .env("PATH", "")
                .current_dir(moved_parent.path())
                .output()
                .unwrap()
        };
        let output = run();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let summary: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(summary["scope"], "selected_files");
        assert_eq!(summary["verified_files"], 2);
        assert_eq!(
            summary["complete_payload_bytes"],
            record.inventory.iter().map(|f| f.bytes).sum::<u64>()
        );
        assert_eq!(
            summary["native_filename_suffixes"][".safetensors"]["files"],
            1
        );
        let payload = moved.join("weights/model.safetensors");
        let mut changed = fs::read(&payload).unwrap();
        changed[299_999] ^= 1;
        fs::write(payload, changed).unwrap();
        let modified = run();
        assert!(!modified.status.success());
        assert!(String::from_utf8_lossy(&modified.stderr).contains("SHA-256 mismatch"));
        // Independent schema validation also rejects a false complete-model claim.
        change_record(&moved, |value| value["scope"] = "complete_model".into());
        let scope = run();
        assert!(!scope.status.success());
        assert!(String::from_utf8_lossy(&scope.stderr).contains("record schema or scope"));
    }
}
