use crate::{
    native::{self, NativeMetadata},
    tree, *,
};
use fs4::fs_std::FileExt;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File, OpenOptions},
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

const MAX_RESPONSE_BYTES: usize = 64 * 1024;
const ROOT_README: &str = "# Portable native biological packages\n\nBrowse artifacts/<namespace>/<canonical-accession>/snapshots/<source-content-id>/.\nEach snapshot README lists ordinary native analysis files in source/. The complete\nprovider tree is preserved; acquisition.json is a small operational receipt and\nchecksums.sha256 verifies saved source bytes with ordinary sha256sum. Native\nREADME/catalog/records remain authoritative. Copy complete snapshots freely;\nno BioV, original import path, catalog database or network is needed to read them.\n.staging is incomplete work, never a ready snapshot. .locks coordinates writers.\nNative registration does not create analysis indices, download data, maintain\naliases, garbage collect or establish provider freshness. Separately prepared\nviews may live under prepared/; each view README declares its native dependency\nclosure and ordinary-reader instructions.\n";

/// A read-only handle to an existing, trusted store root. Construction and
/// resolution do not create directories. Registration alone writes new bundles.
/// Trusted roots are not a sandbox against hostile concurrent writers.
#[derive(Debug, Clone)]
pub struct NativeStore {
    root: PathBuf,
}

impl NativeStore {
    pub fn new(store_root: impl AsRef<Path>) -> Result<Self, StorageError> {
        Ok(Self {
            root: tree::root(store_root.as_ref())?,
        })
    }

    /// Stream a complete copy, validate native metadata, and publish without
    /// replacement. No originals are moved, deleted, mutated or hardlinked.
    pub fn register(
        &self,
        source_root: impl AsRef<Path>,
        request: RegisterRequest,
    ) -> Result<Registration, StorageError> {
        tree::no_links(&self.root)?;
        let source_root = tree::root(source_root.as_ref())?;
        if source_root.starts_with(&self.root) || self.root.starts_with(&source_root) {
            return Err(tree::invalid(
                "configured source and store roots must be disjoint",
            ));
        }
        let source = source_root.join(tree::relative(&request.source_path)?);
        tree::no_links(&source)?;
        let canonical = native::parse_ref(&request.canonical_ref, true)?;
        let requested = native::parse_ref(&request.requested_ref, false)?;
        if !requested.matches(&canonical) {
            return Err(tree::invalid(
                "requested reference does not match the canonical biological reference",
            ));
        }
        let staging_root = tree::ensure_dirs(&self.root, ".staging")?;
        let stage = tempfile::Builder::new()
            .prefix("register-")
            .tempdir_in(&staging_root)
            .map_err(|e| io("create registration staging", e))?;
        let staged_source = stage.path().join("source");
        fs::create_dir(&staged_source).map_err(|e| io("create staged source", e))?;
        let copied = tree::inventory(&source, Some(&staged_source))?;
        // A second complete inventory catches ordinary in-flight edits, additions
        // and removals. Trusted permissions still exclude adversarial TOCTOU.
        if copied != tree::inventory(&source, None)? {
            return Err(tree::package(
                "source changed while registration was in progress",
            ));
        }
        let metadata = native::validate(&staged_source, &canonical, &request.declaration, &copied)?;
        let digest = tree::content_digest(&copied)?;
        let snapshot_id = format!("sha256-{digest}");
        let snapshot_path = format!(
            "artifacts/{}/{}/snapshots/{snapshot_id}",
            canonical.namespace, canonical.accession
        );
        let declaration = normalized_declaration(&request.declaration, &metadata);
        let mut receipt = Receipt {
            schema_version: 1,
            snapshot_id,
            source_content_sha256: digest,
            readme_sha256: String::new(),
            requested_ref: request.requested_ref.clone(),
            canonical_ref: canonical.compact(),
            scope: metadata.scope,
            declaration,
            native_metadata: metadata.metadata,
            representations: metadata.representations,
            unavailable_representations: metadata.unavailable,
            inventory: copied,
            validation: metadata.validation,
            origin: "local_copy".into(),
            registered_at_unix_seconds: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map_err(|_| tree::invalid("system clock is before Unix epoch"))?
                .as_secs(),
            acquired_at: None,
            source_url: None,
            acquisition_tool: None,
        };
        let readme = snapshot_readme(&receipt);
        receipt.readme_sha256 = format!("{:x}", Sha256::digest(readme.as_bytes()));
        let json = serde_json::to_vec_pretty(&receipt)
            .map_err(|_| tree::package("receipt serialization failed"))?;
        tree::write_new(&stage.path().join("acquisition.json"), &json)?;
        tree::write_new(
            &stage.path().join("checksums.sha256"),
            tree::checksum_text(&receipt.inventory).as_bytes(),
        )?;
        tree::write_new(&stage.path().join("README.md"), readme.as_bytes())?;
        sync_tree_directories(&staged_source, &receipt.inventory)?;
        sync_dir(stage.path())?;

        let locks = tree::ensure_dirs(&self.root, ".locks")?;
        let lock_path = locks.join(format!(
            "{}--{}.lock",
            canonical.namespace, canonical.accession
        ));
        match fs::symlink_metadata(&lock_path) {
            Ok(meta) if meta.is_file() && !meta.file_type().is_symlink() => (),
            Ok(_) => {
                return Err(tree::invalid(
                    "registration lock must be a plain regular file",
                ))
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => (),
            Err(e) => return Err(io("inspect registration lock", e)),
        }
        let lock = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(&lock_path)
            .map_err(|e| io("open registration lock", e))?;
        tree::no_links(&lock_path)?;
        if !lock
            .metadata()
            .map_err(|e| io("inspect registration lock", e))?
            .is_file()
        {
            return Err(tree::invalid("registration lock must be a regular file"));
        }
        lock.lock_exclusive()
            .map_err(|e| io("lock snapshot publication", e))?;
        let parent_relative = format!(
            "artifacts/{}/{}/snapshots",
            canonical.namespace, canonical.accession
        );
        let parent = tree::ensure_dirs(&self.root, &parent_relative)?;
        self.ensure_root_readme()?;
        let target = self.root.join(&snapshot_path);
        if target
            .try_exists()
            .map_err(|e| io("inspect target snapshot", e))?
        {
            let winner = self.verify_snapshot(&snapshot_path)?;
            if winner.scope != receipt.scope
                || winner.declaration != receipt.declaration
                || winner.representations != receipt.representations
            {
                return Err(StorageError::DeclarationConflict);
            }
            return bounded(Registration {
                reused: true,
                requested_ref: request.requested_ref,
                snapshot: summary(&snapshot_path, &winner)?,
            });
        }
        // Detect dangling symlinks too; existence alone follows links.
        match fs::symlink_metadata(&target) {
            Ok(_) => {
                return Err(StorageError::Corrupt(
                    "target snapshot already exists".into(),
                ))
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => (),
            Err(e) => return Err(io("inspect publication target", e)),
        }
        // Validate response bounds BEFORE publishing a success that cannot be returned.
        let response = bounded(Registration {
            reused: false,
            requested_ref: request.requested_ref,
            snapshot: summary(&snapshot_path, &receipt)?,
        })?;
        publish_no_replace(stage.path(), &target)?;
        sync_dir(&parent)?;
        sync_dir(&staging_root)?;
        Ok(response)
    }

    /// Resolve only local snapshots. Every matching selected snapshot is fully
    /// rehashed and native validation repeated; there is no network or latest
    /// heuristic. An absent representation is not an absent biological feature.
    pub fn resolve(&self, request: ResolveRequest) -> Result<Resolution, StorageError> {
        tree::no_links(&self.root)?;
        let reference = native::parse_ref(&request.reference, false)?;
        native::representation_name(&request.representation)?;
        let scope = request
            .scope
            .as_deref()
            .unwrap_or(if reference.namespace == "pdb" {
                "entry"
            } else {
                "assembly"
            });
        native::validate_scope(reference.namespace, scope)?;
        if let Some(id) = &request.snapshot_id {
            validate_snapshot_id(id)?;
        }
        let namespace_relative = format!("artifacts/{}", reference.namespace);
        let Some(namespace_dir) = tree::existing_dir(&self.root, &namespace_relative)? else {
            return Ok(Resolution::Miss);
        };
        let mut accessions = Vec::new();
        if reference.version.is_some() || reference.namespace == "pdb" {
            // Canonical references address readable directories directly; unrelated
            // accessions do not consume the bounded discovery scan.
            accessions.push(reference.accession.clone());
        } else {
            for (scanned, entry) in fs::read_dir(&namespace_dir)
                .map_err(|e| io("scan local accessions", e))?
                .enumerate()
            {
                if scanned >= MAX_SCAN_CANDIDATES {
                    return Err(tree::limit("unversioned namespace scan exceeds 1024 entries; specify a biological version"));
                }
                let entry = entry.map_err(|e| io("read accession directory", e))?;
                let Some(name) = entry.file_name().to_str().map(str::to_owned) else {
                    continue;
                };
                let Ok(canonical) =
                    native::parse_ref(&format!("{}:{name}", reference.namespace), true)
                else {
                    continue;
                };
                if reference.matches(&canonical) {
                    accessions.push(name);
                }
            }
        }
        let mut paths = Vec::new();
        let mut scanned_snapshots = 0usize;
        for name in accessions {
            let snapshots_relative = format!("{namespace_relative}/{name}/snapshots");
            let Some(snapshots) = tree::existing_dir(&self.root, &snapshots_relative)? else {
                continue;
            };
            if let Some(pin) = &request.snapshot_id {
                // Exact pins avoid scanning ignored entries or other snapshots.
                match fs::symlink_metadata(snapshots.join(pin)) {
                    Ok(_) => paths.push(format!("{snapshots_relative}/{pin}")),
                    Err(e) if e.kind() == std::io::ErrorKind::NotFound => (),
                    Err(e) => return Err(io("inspect pinned snapshot", e)),
                }
                continue;
            }
            for snapshot in fs::read_dir(&snapshots).map_err(|e| io("scan snapshots", e))? {
                scanned_snapshots += 1;
                if scanned_snapshots > MAX_SCAN_CANDIDATES {
                    return Err(tree::limit(
                        "resolution scan exceeds 1024 directory entries; pin a snapshot",
                    ));
                }
                let snapshot = snapshot.map_err(|e| io("read snapshot directory", e))?;
                let Some(id) = snapshot.file_name().to_str().map(str::to_owned) else {
                    return Err(tree::invalid("snapshot names must be UTF-8"));
                };
                if id.starts_with("sha256-") {
                    paths.push(format!("{snapshots_relative}/{id}"));
                }
            }
        }
        paths.sort();
        let mut candidates = Vec::new();
        let mut selected = None;
        for path in paths {
            let id = path.rsplit('/').next().unwrap_or("").to_owned();
            let receipt = match self.verify_snapshot(&path) {
                Ok(value) => value,
                Err(error) => {
                    return bounded(Resolution::Corrupt {
                        snapshot_id: id,
                        detail: error.to_string(),
                    })
                }
            };
            if receipt.scope != scope {
                continue;
            }
            let candidate = summary(&path, &receipt)?;
            if candidates.len() >= MAX_AMBIGUOUS_CANDIDATES {
                return Err(tree::limit(
                    "ambiguity exceeds 64 snapshots; pin a source content ID",
                ));
            }
            candidates.push(candidate);
            if candidates.len() == 1 {
                selected = Some((path, receipt));
            }
        }
        if candidates.is_empty() {
            return Ok(Resolution::Miss);
        }
        if candidates.len() > 1 {
            return bounded(Resolution::Ambiguous { candidates });
        }
        let (path, receipt) =
            selected.ok_or_else(|| StorageError::Corrupt("missing selected candidate".into()))?;
        let snapshot = candidates.remove(0);
        let Some(files) = receipt.representations.get(&request.representation) else {
            return bounded(Resolution::Unavailable {
                snapshot,
                available_representations: receipt.representations.keys().cloned().collect(),
            });
        };
        let paths = files
            .iter()
            .map(|relative| {
                let relative_path = format!("{path}/source/{relative}");
                let execution_host_path = self
                    .root
                    .join(&relative_path)
                    .to_str()
                    .ok_or_else(|| tree::invalid("execution host path is not UTF-8"))?
                    .to_owned();
                Ok(ResolvedFile {
                    relative_path,
                    execution_host_path,
                })
            })
            .collect::<Result<Vec<_>, StorageError>>()?;
        bounded(Resolution::Ready { snapshot, paths })
    }

    fn ensure_root_readme(&self) -> Result<(), StorageError> {
        let readme = self.root.join("README.md");
        match fs::symlink_metadata(&readme) {
            Ok(meta) if meta.is_file() && !meta.file_type().is_symlink() => (), // Preserve a user's existing root documentation.
            Ok(_) => return Err(tree::invalid("store README is not a regular file")),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                match tree::write_new(&readme, ROOT_README.as_bytes()) {
                    Ok(()) => (),
                    Err(StorageError::Io {
                        kind: std::io::ErrorKind::AlreadyExists,
                        ..
                    }) => (),
                    Err(error) => return Err(error),
                }
            }
            Err(e) => return Err(io("inspect root README", e)),
        }
        tree::no_links(&readme)?;
        if !fs::metadata(&readme)
            .map_err(|e| io("inspect root README", e))?
            .is_file()
        {
            return Err(tree::invalid("store README is not a regular file"));
        }
        sync_dir(&self.root)
    }

    fn verify_snapshot(&self, path: &str) -> Result<Receipt, StorageError> {
        let snapshot = self.root.join(tree::relative(path)?);
        tree::no_links(&snapshot)?;
        let receipt: Receipt =
            serde_json::from_slice(&tree::read_bounded(&snapshot.join("acquisition.json"))?)
                .map_err(|_| StorageError::Corrupt("malformed acquisition receipt".into()))?;
        if receipt.schema_version != 1
            || receipt.origin != "local_copy"
            || receipt.acquired_at.is_some()
            || receipt.source_url.is_some()
            || receipt.acquisition_tool.is_some()
        {
            return Err(StorageError::Corrupt(
                "unsupported receipt schema or acquisition claims".into(),
            ));
        }
        validate_snapshot_id(&receipt.snapshot_id)?;
        let canonical = native::parse_ref(&receipt.canonical_ref, true)?;
        let requested = native::parse_ref(&receipt.requested_ref, false)?;
        if canonical.compact() != receipt.canonical_ref
            || !requested.matches(&canonical)
            || path
                != format!(
                    "artifacts/{}/{}/snapshots/{}",
                    canonical.namespace, canonical.accession, receipt.snapshot_id
                )
        {
            return Err(StorageError::Corrupt(
                "receipt reference/path mismatch".into(),
            ));
        }
        let digest = tree::content_digest(&receipt.inventory)?;
        if digest != receipt.source_content_sha256
            || receipt.snapshot_id != format!("sha256-{digest}")
        {
            return Err(StorageError::Corrupt(
                "source inventory identity mismatch".into(),
            ));
        }
        let actual = tree::inventory(&snapshot.join("source"), None)?;
        if actual != receipt.inventory {
            return Err(StorageError::Corrupt(
                "native source tree bytes or paths changed".into(),
            ));
        }
        let native = native::validate(
            &snapshot.join("source"),
            &canonical,
            &receipt.declaration,
            &actual,
        )?;
        if native.scope != receipt.scope
            || native.metadata != receipt.native_metadata
            || native.representations != receipt.representations
            || native.unavailable != receipt.unavailable_representations
            || native.validation.method != receipt.validation.method
            || normalized_declaration(&receipt.declaration, &native) != receipt.declaration
        {
            return Err(StorageError::Corrupt(
                "receipt disagrees with native metadata or declaration".into(),
            ));
        }
        if tree::read_bounded(&snapshot.join("checksums.sha256"))?
            != tree::checksum_text(&receipt.inventory).as_bytes()
        {
            return Err(StorageError::Corrupt(
                "checksum listing disagrees with native inventory".into(),
            ));
        }
        let readme = tree::read_bounded(&snapshot.join("README.md"))?;
        if format!("{:x}", Sha256::digest(&readme)) != receipt.readme_sha256 {
            return Err(StorageError::Corrupt(
                "snapshot entry-point checksum mismatch".into(),
            ));
        }
        let mut wrappers = 0;
        for member in fs::read_dir(&snapshot).map_err(|e| io("inspect snapshot wrapper", e))? {
            let member = member.map_err(|e| io("read snapshot wrapper", e))?;
            wrappers += 1;
            if !matches!(
                member.file_name().to_str(),
                Some("source" | "acquisition.json" | "checksums.sha256" | "README.md")
            ) {
                return Err(StorageError::Corrupt(
                    "unexpected snapshot wrapper member".into(),
                ));
            }
        }
        if wrappers != 4 {
            return Err(StorageError::Corrupt("incomplete snapshot wrapper".into()));
        }
        Ok(receipt)
    }
}

fn normalized_declaration(
    declaration: &NativeDeclaration,
    metadata: &NativeMetadata,
) -> NativeDeclaration {
    match declaration {
        NativeDeclaration::Refseq => NativeDeclaration::Refseq,
        NativeDeclaration::Pdb { .. } => NativeDeclaration::Pdb {
            scope: metadata.scope.clone(),
            representations: metadata.representations.clone(),
        },
    }
}
fn validate_snapshot_id(id: &str) -> Result<(), StorageError> {
    let digest = id.strip_prefix("sha256-").ok_or_else(|| {
        tree::invalid("snapshot ID must be sha256- plus 64 lowercase hexadecimal characters")
    })?;
    tree::decode_hash(digest).map_err(|_| tree::invalid("invalid snapshot ID"))?;
    Ok(())
}
fn summary(path: &str, receipt: &Receipt) -> Result<SnapshotSummary, StorageError> {
    let total_bytes = receipt
        .inventory
        .iter()
        .try_fold(0u64, |sum, e| sum.checked_add(e.bytes))
        .ok_or_else(|| tree::limit("total byte count overflow"))?;
    Ok(SnapshotSummary {
        snapshot_id: receipt.snapshot_id.clone(),
        snapshot_path: path.into(),
        receipt_path: format!("{path}/acquisition.json"),
        source_content_sha256: receipt.source_content_sha256.clone(),
        requested_ref: receipt.requested_ref.clone(),
        canonical_ref: receipt.canonical_ref.clone(),
        scope: receipt.scope.clone(),
        file_count: receipt
            .inventory
            .iter()
            .filter(|e| e.kind == EntryKind::File)
            .count(),
        total_bytes,
    })
}
fn bounded<T: Serialize>(value: T) -> Result<T, StorageError> {
    let bytes =
        serde_json::to_vec(&value).map_err(|_| tree::limit("cannot serialize storage response"))?;
    if bytes.len() > MAX_RESPONSE_BYTES {
        return Err(tree::limit(
            "storage response exceeds 64 KiB; use a narrower selection",
        ));
    }
    Ok(value)
}
fn markdown_code(value: &str) -> String {
    let longest = value.split(|c| c != '`').map(str::len).max().unwrap_or(0);
    let delimiter = "`".repeat(longest + 1);
    // CommonMark removes one padding space on each side, preserving backticks
    // and any meaningful internal spaces in these nonempty relative paths.
    format!("{delimiter} {value} {delimiter}")
}
fn snapshot_readme(receipt: &Receipt) -> String {
    let mut text = format!("# Native package: {}\n\nScope: {}. Origin: copied local directory; upstream acquisition time, URL and\ntool are unknown. Registration time in acquisition.json is not acquisition time.\nSnapshot: {}. This identifies native source bytes and paths, not a biological\nrelease/version. Native metadata in source/ remains authoritative.\n\n## Analysis entry files\n\n", receipt.canonical_ref, receipt.scope, receipt.snapshot_id);
    for (name, files) in &receipt.representations {
        for file in files {
            text.push_str(&format!(
                "- `{name}`: {}\n",
                markdown_code(&format!("source/{file}"))
            ));
        }
    }
    if !receipt.unavailable_representations.is_empty() {
        text.push_str(&format!("\nNot present in this native package: {}. This is local representation\navailability, not evidence of biological absence.\n", receipt.unavailable_representations.join(", ")));
    }
    text.push_str("\n## Native documentation and integrity\n\n");
    if receipt.native_metadata.is_empty() {
        text.push_str("Identity, scope and representation mapping were supplied by the caller.\nNo provider identity, structure parsing, assembly identity or scientific QC was verified.\n");
    }
    for path in &receipt.native_metadata {
        text.push_str(&format!("- {}\n", markdown_code(&format!("source/{path}"))));
    }
    text.push_str("\nCopy this entire snapshot directory anywhere. From its root, run:\n\n```sh\nsha256sum -c checksums.sha256\n```\n\nThe portable JSON receipt contains sorted source-relative paths, sizes and SHA-256\nhashes; directory entries preserve empty directories. Hashes describe saved file\nbytes (compressed bytes for retained .gz), not decoded biological records. Checksums\nestablish consistency with supplied records, not producer authenticity or scientific QC.\nThere are no derived indices here. Full-source hashing is performed on explicit\nregistration and resolution; no large-data performance claim is implied.\n\n## Ordinary-reader example\n\n");
    if let Some(paths) = receipt.representations.get("genome_fasta") {
        let paths: Vec<_> = paths.iter().map(|path| format!("source/{path}")).collect();
        let paths = serde_json::to_string(&paths).unwrap_or_default();
        text.push_str(&format!("Use any ordinary FASTA reader; for example, Python's standard library counts\ncomplete records and bases across every genome FASTA member, without BioV:\n\n```python\nfrom pathlib import Path\npaths = {paths}\nrecords = bases = 0\nfor p in map(Path, paths):\n    with p.open() as handle:\n        for line in handle:\n            if line.startswith('>'):\n                records += 1\n            else:\n                bases += len(line.strip())\nprint(records, bases)\n```\n\nThe genome and listed annotation files are from this same native package. Assembly\naccession versions do not freeze future annotation/package bytes.\n"));
    } else if let Some(path) = receipt
        .representations
        .values()
        .flatten()
        .find(|p| p.ends_with(".cif") || p.ends_with(".cif.gz"))
    {
        let path = serde_json::to_string(&format!("source/{path}")).unwrap_or_default();
        text.push_str(&format!("For one listed caller-declared coordinate file, if it is mmCIF, an ordinary\nBiopython reader can read every record in that file (install Biopython separately). This example does not establish identity:\n\n```python\nimport gzip\nfrom pathlib import Path\nfrom Bio.PDB import MMCIFParser\np = Path({path})\nopener = gzip.open if p.suffix == '.gz' else open\nwith opener(p, 'rt') as handle:\n    structure = MMCIFParser(QUIET=True).get_structure('declared', handle)\nprint(sum(1 for _ in structure.get_residues()), sum(1 for _ in structure.get_atoms()))\n```\n"));
    } else {
        text.push_str("Open the entry files above with ordinary readers for their native formats; the\ncaller-declared mapping does not validate formats or biological identity.\n");
    }
    text
}

#[cfg(any(target_os = "linux", target_os = "android", target_vendor = "apple"))]
fn publish_no_replace(source: &Path, target: &Path) -> Result<(), StorageError> {
    rustix::fs::renameat_with(
        rustix::fs::CWD,
        source,
        rustix::fs::CWD,
        target,
        rustix::fs::RenameFlags::NOREPLACE,
    )
    .map_err(|e| {
        if e == rustix::io::Errno::NOSYS
            || e == rustix::io::Errno::NOTSUP
            || e == rustix::io::Errno::INVAL
        {
            StorageError::UnsupportedPublication
        } else {
            io(
                "atomic no-replace snapshot publication",
                std::io::Error::from(e),
            )
        }
    })
}
#[cfg(windows)]
fn publish_no_replace(source: &Path, target: &Path) -> Result<(), StorageError> {
    fs::rename(source, target).map_err(|e| io("atomic no-replace snapshot publication", e))
}
#[cfg(not(any(
    target_os = "linux",
    target_os = "android",
    target_vendor = "apple",
    windows
)))]
fn publish_no_replace(_: &Path, _: &Path) -> Result<(), StorageError> {
    Err(StorageError::UnsupportedPublication)
}
#[cfg(unix)]
fn sync_dir(path: &Path) -> Result<(), StorageError> {
    File::open(path)
        .and_then(|f| f.sync_all())
        .map_err(|e| io("sync directory", e))
}
#[cfg(not(unix))]
fn sync_dir(_: &Path) -> Result<(), StorageError> {
    Ok(())
}
fn sync_tree_directories(source: &Path, inventory: &[InventoryEntry]) -> Result<(), StorageError> {
    for entry in inventory
        .iter()
        .rev()
        .filter(|e| e.kind == EntryKind::Directory)
    {
        sync_dir(&source.join(&entry.path))?;
    }
    sync_dir(source)
}
