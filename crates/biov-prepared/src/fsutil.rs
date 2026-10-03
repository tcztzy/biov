//! Local publication safety, under the same trusted-root contract as biov-storage.
use crate::{corrupt, invalid, io, limit, OutputIdentity, PreparedError, IO_BUFFER_BYTES};
use fs4::fs_std::FileExt;
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Component, Path, PathBuf},
};

pub(crate) fn relative(value: &str) -> Result<PathBuf, PreparedError> {
    if value.is_empty()
        || value.len() > 4096
        || value.contains(['\\', ':'])
        || value.chars().any(char::is_control)
        || value.starts_with('/')
        || value
            .split('/')
            .any(|v| v.is_empty() || v == "." || v == ".." || v.ends_with([' ', '.']))
    {
        return Err(invalid("expected a bounded normal portable relative path"));
    }
    for part in value.split('/') {
        let stem = part.split('.').next().unwrap_or("").to_ascii_uppercase();
        if matches!(stem.as_str(), "CON" | "PRN" | "AUX" | "NUL")
            || (stem.len() == 4
                && (stem.starts_with("COM") || stem.starts_with("LPT"))
                && matches!(stem.as_bytes()[3], b'1'..=b'9'))
            || part.contains(['<', '>', '"', '|', '?', '*'])
        {
            return Err(invalid("source path is not portable"));
        }
    }
    let path = Path::new(value);
    if !path
        .components()
        .all(|part| matches!(part, Component::Normal(_)))
    {
        return Err(invalid("source path contains non-normal components"));
    }
    Ok(path.to_owned())
}

pub(crate) fn no_links(path: &Path) -> Result<(), PreparedError> {
    let mut current = PathBuf::new();
    for part in path.components() {
        if matches!(part, Component::ParentDir) {
            return Err(invalid("parent path components are not allowed"));
        }
        current.push(part);
        let metadata = fs::symlink_metadata(&current).map_err(|e| io("inspect path", e))?;
        if metadata.file_type().is_symlink() {
            return Err(invalid("symbolic links are not allowed"));
        }
    }
    Ok(())
}

pub(crate) fn root(path: &Path) -> Result<PathBuf, PreparedError> {
    if path.to_str().is_none() {
        return Err(invalid("store root must be UTF-8"));
    }
    let absolute = if path.is_absolute() {
        path.to_owned()
    } else {
        std::env::current_dir()
            .map_err(|e| io("read working directory", e))?
            .join(path)
    };
    no_links(&absolute)?;
    let root = fs::canonicalize(absolute).map_err(|e| io("canonicalize root", e))?;
    if !root.is_dir() {
        return Err(invalid("store root must be an existing directory"));
    }
    Ok(root)
}

pub(crate) fn ensure_dir(root: &Path, name: &str) -> Result<PathBuf, PreparedError> {
    let mut path = root.to_owned();
    for part in relative(name)?.components() {
        path.push(part);
        match fs::create_dir(&path) {
            Ok(()) => sync_dir(path.parent().ok_or_else(|| invalid("missing parent"))?)?,
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => (),
            Err(e) => return Err(io("create prepared directory", e)),
        }
        let metadata =
            fs::symlink_metadata(&path).map_err(|e| io("inspect prepared directory", e))?;
        if !metadata.is_dir() || metadata.file_type().is_symlink() {
            return Err(invalid("store component is not a plain directory"));
        }
    }
    Ok(path)
}

pub(crate) fn open_regular(path: &Path) -> Result<File, PreparedError> {
    no_links(path)?;
    if !fs::metadata(path)
        .map_err(|e| io("inspect file", e))?
        .is_file()
    {
        return Err(invalid("expected a plain regular file"));
    }
    let file = File::open(path).map_err(|e| io("open regular file", e))?;
    if !file
        .metadata()
        .map_err(|e| io("inspect opened file", e))?
        .is_file()
    {
        return Err(invalid("opened file is not regular"));
    }
    Ok(file)
}

pub(crate) fn read_bounded(path: &Path, max: usize) -> Result<Vec<u8>, PreparedError> {
    let file = open_regular(path)?;
    if file
        .metadata()
        .map_err(|e| io("inspect metadata", e))?
        .len()
        > max as u64
    {
        return Err(limit("metadata document exceeds its byte limit"));
    }
    let mut bytes = Vec::new();
    file.take((max + 1) as u64)
        .read_to_end(&mut bytes)
        .map_err(|e| io("read metadata", e))?;
    if bytes.len() > max {
        return Err(limit("metadata document exceeds its byte limit"));
    }
    Ok(bytes)
}

pub(crate) fn create_file(path: &Path) -> Result<File, PreparedError> {
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| io("create prepared file", e))
}
pub(crate) fn write_new(path: &Path, bytes: &[u8]) -> Result<(), PreparedError> {
    let mut file = create_file(path)?;
    file.write_all(bytes)
        .map_err(|e| io("write prepared file", e))?;
    file.sync_all().map_err(|e| io("sync prepared file", e))
}

pub(crate) fn hash_file(
    path: &Path,
    name: &str,
    max: usize,
) -> Result<OutputIdentity, PreparedError> {
    let mut file = open_regular(path)?;
    if file
        .metadata()
        .map_err(|e| io("inspect output size", e))?
        .len()
        > max as u64
    {
        return Err(limit("output exceeds metadata byte limit"));
    }
    let mut digest = Sha256::new();
    let mut count = 0u64;
    let mut buf = [0; IO_BUFFER_BYTES];
    loop {
        let n = file
            .read(&mut buf)
            .map_err(|e| io("hash prepared output", e))?;
        if n == 0 {
            break;
        }
        count += n as u64;
        if count > max as u64 {
            return Err(limit("output exceeds metadata byte limit"));
        }
        digest.update(&buf[..n]);
    }
    Ok(OutputIdentity {
        path: name.into(),
        bytes: count,
        sha256: format!("{:x}", digest.finalize()),
    })
}

pub(crate) fn lock(root: &Path, recipe_id: &str) -> Result<File, PreparedError> {
    let path = ensure_dir(root, ".locks")?.join(format!("prepared-{recipe_id}.lock"));
    match fs::symlink_metadata(&path) {
        Ok(meta) if meta.is_file() && !meta.file_type().is_symlink() => (),
        Ok(_) => return Err(invalid("prepared lock must be a plain regular file")),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => (),
        Err(e) => return Err(io("inspect prepared lock", e)),
    }
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(&path)
        .map_err(|e| io("open prepared lock", e))?;
    no_links(&path)?;
    if !file
        .metadata()
        .map_err(|e| io("inspect prepared lock", e))?
        .is_file()
    {
        return Err(invalid("prepared lock must be regular"));
    }
    file.lock_exclusive()
        .map_err(|e| io("lock prepared publication", e))?;
    Ok(file)
}

pub(crate) fn verify_members(path: &Path) -> Result<(), PreparedError> {
    no_links(path)?;
    if !path.is_dir() {
        return Err(corrupt("prepared target is not a directory"));
    }
    let mut count = 0;
    for member in fs::read_dir(path).map_err(|e| io("read prepared directory", e))? {
        let member = member.map_err(|e| io("read prepared member", e))?;
        count += 1;
        if count > 4
            || !matches!(
                member.file_name().to_str(),
                Some("sequences.fai" | "sequences.tsv" | "README.md" | "provenance.json")
            )
        {
            return Err(corrupt("unexpected prepared directory member"));
        }
        let metadata = member
            .file_type()
            .map_err(|e| io("inspect prepared member", e))?;
        if !metadata.is_file() {
            return Err(corrupt("prepared members must be plain regular files"));
        }
    }
    if count != 4 {
        return Err(corrupt("incomplete prepared directory"));
    }
    Ok(())
}

#[cfg(any(target_os = "linux", target_os = "android", target_vendor = "apple"))]
pub(crate) fn publish(source: &Path, target: &Path) -> Result<(), PreparedError> {
    rustix::fs::renameat_with(
        rustix::fs::CWD,
        source,
        rustix::fs::CWD,
        target,
        rustix::fs::RenameFlags::NOREPLACE,
    )
    .map_err(|e| {
        if matches!(
            e,
            rustix::io::Errno::NOSYS | rustix::io::Errno::NOTSUP | rustix::io::Errno::INVAL
        ) {
            PreparedError::UnsupportedPublication
        } else {
            io(
                "atomic no-replace prepared publication",
                std::io::Error::from(e),
            )
        }
    })
}
#[cfg(windows)]
pub(crate) fn publish(source: &Path, target: &Path) -> Result<(), PreparedError> {
    fs::rename(source, target).map_err(|e| io("atomic no-replace prepared publication", e))
}
#[cfg(not(any(
    target_os = "linux",
    target_os = "android",
    target_vendor = "apple",
    windows
)))]
pub(crate) fn publish(_: &Path, _: &Path) -> Result<(), PreparedError> {
    Err(PreparedError::UnsupportedPublication)
}
#[cfg(unix)]
pub(crate) fn sync_dir(path: &Path) -> Result<(), PreparedError> {
    File::open(path)
        .and_then(|f| f.sync_all())
        .map_err(|e| io("sync prepared directory", e))
}
#[cfg(not(unix))]
pub(crate) fn sync_dir(_: &Path) -> Result<(), PreparedError> {
    Ok(())
}
