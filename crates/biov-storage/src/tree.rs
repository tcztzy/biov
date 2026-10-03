use crate::{
    io, EntryKind, InventoryEntry, StorageError, MAX_ENTRIES, MAX_METADATA_BYTES, MAX_PATH_BYTES,
};
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Component, Path, PathBuf},
};

pub(crate) fn invalid(message: &str) -> StorageError {
    StorageError::InvalidInput(message.into())
}
pub(crate) fn package(message: &str) -> StorageError {
    StorageError::InvalidPackage(message.into())
}
pub(crate) fn limit(message: &str) -> StorageError {
    StorageError::Limit(message.into())
}

pub(crate) fn relative(value: &str) -> Result<PathBuf, StorageError> {
    if value.is_empty()
        || value.len() > MAX_PATH_BYTES
        || value.contains('\\')
        || value.contains(':')
        || value.chars().any(char::is_control)
        || value.starts_with('/')
        || value
            .split('/')
            .any(|v| v.is_empty() || v == "." || v == ".." || v.ends_with([' ', '.']))
    {
        return Err(invalid(
            "paths must be bounded, portable, normal UTF-8 relative paths",
        ));
    }
    let path = Path::new(value);
    if !path
        .components()
        .all(|part| matches!(part, Component::Normal(_)))
    {
        return Err(invalid("path contains a non-normal component"));
    }
    // Reject Windows device names on every platform so moved bundles stay usable.
    for part in value.split('/') {
        let stem = part.split('.').next().unwrap_or("").to_ascii_uppercase();
        if matches!(stem.as_str(), "CON" | "PRN" | "AUX" | "NUL")
            || (stem.len() == 4
                && (stem.starts_with("COM") || stem.starts_with("LPT"))
                && matches!(stem.as_bytes()[3], b'1'..=b'9'))
            || part.contains(['<', '>', '"', '|', '?', '*'])
        {
            return Err(invalid("path is not portable across supported filesystems"));
        }
    }
    Ok(path.to_owned())
}

pub(crate) fn no_links(path: &Path) -> Result<(), StorageError> {
    let mut current = PathBuf::new();
    for part in path.components() {
        if matches!(part, Component::ParentDir) {
            return Err(invalid("parent path components are not allowed"));
        }
        current.push(part);
        let meta = fs::symlink_metadata(&current).map_err(|e| io("inspect path", e))?;
        if meta.file_type().is_symlink() {
            return Err(invalid("symbolic links are not allowed"));
        }
    }
    Ok(())
}

pub(crate) fn root(path: &Path) -> Result<PathBuf, StorageError> {
    if path.to_str().is_none() {
        return Err(invalid("root paths must be UTF-8"));
    }
    let absolute = if path.is_absolute() {
        path.to_owned()
    } else {
        std::env::current_dir()
            .map_err(|e| io("read working directory", e))?
            .join(path)
    };
    no_links(&absolute)?;
    let canonical = fs::canonicalize(absolute).map_err(|e| io("canonicalize root", e))?;
    if !canonical.is_dir() {
        return Err(invalid("configured roots must be existing directories"));
    }
    Ok(canonical)
}

pub(crate) fn ensure_dirs(root: &Path, rel: &str) -> Result<PathBuf, StorageError> {
    let rel = relative(rel)?;
    let mut path = root.to_owned();
    for part in rel.components() {
        path.push(part);
        match fs::create_dir(&path) {
            Ok(()) => {
                sync_directory(
                    path.parent()
                        .ok_or_else(|| invalid("directory has no parent"))?,
                )?;
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
            Err(e) => return Err(io("create store directory", e)),
        }
        let meta = fs::symlink_metadata(&path).map_err(|e| io("inspect store directory", e))?;
        if !meta.is_dir() || meta.file_type().is_symlink() {
            return Err(invalid("store component is not a plain directory"));
        }
    }
    Ok(path)
}

pub(crate) fn read_bounded(path: &Path) -> Result<Vec<u8>, StorageError> {
    no_links(path)?;
    let meta = fs::metadata(path).map_err(|e| io("inspect metadata file", e))?;
    if !meta.is_file() {
        return Err(package("metadata must be a regular file"));
    }
    if meta.len() > MAX_METADATA_BYTES as u64 {
        return Err(limit("metadata exceeds 16 MiB"));
    }
    let mut contents = Vec::new();
    File::open(path)
        .map_err(|e| io("open metadata file", e))?
        .take((MAX_METADATA_BYTES + 1) as u64)
        .read_to_end(&mut contents)
        .map_err(|e| io("read metadata file", e))?;
    if contents.len() > MAX_METADATA_BYTES {
        return Err(limit("metadata exceeds 16 MiB"));
    }
    Ok(contents)
}

pub(crate) fn write_new(path: &Path, contents: &[u8]) -> Result<(), StorageError> {
    if contents.len() > MAX_METADATA_BYTES {
        return Err(limit("generated metadata exceeds 16 MiB"));
    }
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| io("create metadata", e))?;
    file.write_all(contents)
        .map_err(|e| io("write metadata", e))?;
    file.sync_all().map_err(|e| io("sync metadata", e))
}

/// Hash all file bytes, optionally copying into a new tree. Never follows links.
pub(crate) fn inventory(
    source: &Path,
    destination: Option<&Path>,
) -> Result<Vec<InventoryEntry>, StorageError> {
    no_links(source)?;
    if !source.is_dir() {
        return Err(invalid("source must be a directory"));
    }
    let mut entries = Vec::new();
    let mut inventory_budget = 0usize;
    let mut pending = vec![(source.to_owned(), String::new())];
    while let Some((dir, prefix)) = pending.pop() {
        for member in fs::read_dir(dir).map_err(|e| io("enumerate native tree", e))? {
            let member = member.map_err(|e| io("read native tree entry", e))?;
            if entries.len() >= MAX_ENTRIES {
                return Err(limit("native tree exceeds 100000 entries"));
            }
            let name = member
                .file_name()
                .into_string()
                .map_err(|_| invalid("native filenames must be UTF-8"))?;
            let rel = if prefix.is_empty() {
                name
            } else {
                format!("{prefix}/{name}")
            };
            relative(&rel)?;
            // Bound retained inventory metadata before copying another source
            // file. This is deliberately independent of file byte counts.
            inventory_budget = inventory_budget
                .checked_add(rel.len() + 256)
                .ok_or_else(|| limit("inventory metadata size overflow"))?;
            if inventory_budget > MAX_METADATA_BYTES {
                return Err(limit("inventory metadata budget exceeds 16 MiB"));
            }
            let metadata =
                fs::symlink_metadata(member.path()).map_err(|e| io("inspect native entry", e))?;
            if metadata.file_type().is_symlink() {
                return Err(invalid("symbolic links are not allowed in native trees"));
            }
            if metadata.is_dir() {
                if let Some(dest) = destination {
                    fs::create_dir(dest.join(&rel)).map_err(|e| io("copy native directory", e))?;
                }
                entries.push(InventoryEntry {
                    path: rel.clone(),
                    kind: EntryKind::Directory,
                    bytes: 0,
                    sha256: None,
                });
                pending.push((member.path(), rel));
            } else if metadata.is_file() {
                let mut input = File::open(member.path()).map_err(|e| io("open native file", e))?;
                let mut output = destination
                    .map(|d| {
                        OpenOptions::new()
                            .write(true)
                            .create_new(true)
                            .open(d.join(&rel))
                    })
                    .transpose()
                    .map_err(|e| io("create native copy", e))?;
                let mut hash = Sha256::new();
                let mut bytes = 0u64;
                let mut buffer = [0u8; 64 * 1024];
                loop {
                    let count = input
                        .read(&mut buffer)
                        .map_err(|e| io("read native bytes", e))?;
                    if count == 0 {
                        break;
                    }
                    hash.update(&buffer[..count]);
                    bytes = bytes
                        .checked_add(count as u64)
                        .ok_or_else(|| limit("file byte count overflow"))?;
                    if let Some(out) = &mut output {
                        out.write_all(&buffer[..count])
                            .map_err(|e| io("write native bytes", e))?;
                    }
                }
                if bytes != metadata.len() {
                    return Err(package("source changed while reading"));
                }
                if let Some(out) = output {
                    out.sync_all().map_err(|e| io("sync native copy", e))?;
                }
                entries.push(InventoryEntry {
                    path: rel,
                    kind: EntryKind::File,
                    bytes,
                    sha256: Some(format!("{:x}", hash.finalize())),
                });
            } else {
                return Err(invalid("special files are not allowed in native trees"));
            }
        }
    }
    entries.sort_by(|a, b| a.path.as_bytes().cmp(b.path.as_bytes()));
    Ok(entries)
}

pub(crate) fn decode_hash(value: &str) -> Result<[u8; 32], StorageError> {
    if value.len() != 64
        || !value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(package("invalid lowercase SHA-256"));
    }
    let mut bytes = [0; 32];
    for (i, byte) in bytes.iter_mut().enumerate() {
        *byte = u8::from_str_radix(&value[i * 2..i * 2 + 2], 16)
            .map_err(|_| package("invalid SHA-256"))?;
    }
    Ok(bytes)
}

pub(crate) fn content_digest(entries: &[InventoryEntry]) -> Result<String, StorageError> {
    if entries.len() > MAX_ENTRIES {
        return Err(limit("inventory exceeds 100000 entries"));
    }
    let mut hash = Sha256::new();
    hash.update(b"biov-native-tree-v1\0");
    let mut previous: Option<&str> = None;
    for entry in entries {
        relative(&entry.path)?;
        if previous.is_some_and(|p| p >= entry.path.as_str()) {
            return Err(package("inventory paths must be unique and byte-sorted"));
        }
        previous = Some(&entry.path);
        let digest = match entry.kind {
            EntryKind::Directory if entry.bytes == 0 && entry.sha256.is_none() => {
                hash.update(b"D");
                [0; 32]
            }
            EntryKind::File => {
                hash.update(b"F");
                decode_hash(
                    entry
                        .sha256
                        .as_deref()
                        .ok_or_else(|| package("file hash missing"))?,
                )?
            }
            _ => return Err(package("invalid directory inventory record")),
        };
        hash.update((entry.path.len() as u32).to_be_bytes());
        hash.update(entry.path.as_bytes());
        hash.update(entry.bytes.to_be_bytes());
        hash.update(digest);
    }
    Ok(format!("{:x}", hash.finalize()))
}

pub(crate) fn checksum_text(entries: &[InventoryEntry]) -> String {
    let mut text = String::new();
    for entry in entries.iter().filter(|e| e.kind == EntryKind::File) {
        text.push_str(entry.sha256.as_deref().unwrap_or(""));
        text.push_str("  source/");
        text.push_str(&entry.path);
        text.push('\n');
    }
    text
}

/// Read-only existence check that rejects links even in intermediate components.
pub(crate) fn existing_dir(root: &Path, rel: &str) -> Result<Option<PathBuf>, StorageError> {
    let rel = relative(rel)?;
    let mut path = root.to_owned();
    for component in rel.components() {
        path.push(component);
        match fs::symlink_metadata(&path) {
            Ok(meta) if meta.is_dir() && !meta.file_type().is_symlink() => (),
            Ok(_) => return Err(invalid("store component is not a plain directory")),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(e) => return Err(io("inspect store directory", e)),
        }
    }
    Ok(Some(path))
}
#[cfg(unix)]
fn sync_directory(path: &Path) -> Result<(), StorageError> {
    File::open(path)
        .and_then(|f| f.sync_all())
        .map_err(|e| io("sync created directory parent", e))
}
#[cfg(not(unix))]
fn sync_directory(_: &Path) -> Result<(), StorageError> {
    Ok(())
}
