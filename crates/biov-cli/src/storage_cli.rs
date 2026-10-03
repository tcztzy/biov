//! One-shot JSON/file adapters; all storage semantics live in `biov-storage`.
use std::{
    fs::File,
    io::{self, Read, Write},
    path::Path,
};

use biov_storage::{NativeStore, RegisterRequest, ResolveRequest};
use serde::{de::DeserializeOwned, Serialize};

const MAX_REQUEST_BYTES: usize = 64 * 1024;
const MAX_RESPONSE_BYTES: usize = 64 * 1024;

pub(crate) fn register(
    store_root: &Path,
    source_root: &Path,
    request_file: &Path,
) -> Result<(), String> {
    let request: RegisterRequest = read_request(request_file)?;
    let store = NativeStore::new(store_root).map_err(|error| error.to_string())?;
    let result = store
        .register(source_root, request)
        .map_err(|error| error.to_string())?;
    write_result(&result)
}

pub(crate) fn resolve(store_root: &Path, request_file: &Path) -> Result<(), String> {
    let request: ResolveRequest = read_request(request_file)?;
    let store = NativeStore::new(store_root).map_err(|error| error.to_string())?;
    // Miss, ambiguity, unavailable representation and corruption are structured
    // outcomes, not process errors. The caller can inspect `status` directly.
    let result = store.resolve(request).map_err(|error| error.to_string())?;
    write_result(&result)
}

pub(crate) fn read_request<T: DeserializeOwned>(path: &Path) -> Result<T, String> {
    // Reject FIFOs/devices before opening: opening a FIFO can block indefinitely.
    // Explicit request-file symlinks are allowed, unlike stored native paths.
    // Retain the handle check below for ordinary changes during opening. As with
    // the configured roots, this is not a hostile-filesystem race sandbox.
    if !std::fs::metadata(path)
        .map_err(|error| format!("cannot inspect request file: {}", error.kind()))?
        .is_file()
    {
        return Err("request file must be a regular file".into());
    }
    let file =
        File::open(path).map_err(|error| format!("cannot open request file: {}", error.kind()))?;
    if !file
        .metadata()
        .map_err(|error| format!("cannot inspect request file: {}", error.kind()))?
        .is_file()
    {
        return Err("request file must be a regular file".into());
    }
    decode_request(file)
}

fn decode_request<T: DeserializeOwned>(reader: impl Read) -> Result<T, String> {
    let mut bytes = Vec::new();
    // Never read the whole input before checking the limit, including if a file
    // changes between metadata inspection and reading.
    reader
        .take((MAX_REQUEST_BYTES + 1) as u64)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("cannot read request file: {}", error.kind()))?;
    if bytes.len() > MAX_REQUEST_BYTES {
        return Err("request file exceeds the 64 KiB limit".into());
    }
    serde_json::from_slice(&bytes).map_err(|error| format!("invalid request JSON: {error}"))
}

pub(crate) fn write_result(value: &impl Serialize) -> Result<(), String> {
    let mut bytes =
        serde_json::to_vec(value).map_err(|error| format!("cannot serialize result: {error}"))?;
    if bytes.len() > MAX_RESPONSE_BYTES {
        return Err("result exceeds the 64 KiB limit".into());
    }
    bytes.push(b'\n');
    let mut stdout = io::stdout().lock();
    stdout
        .write_all(&bytes)
        .and_then(|_| stdout.flush())
        .map_err(|error| format!("cannot write result: {}", error.kind()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const RESOLVE: &str = r#"{"reference":"pdb:1ABC","representation":"mmcif"}"#;

    #[test]
    fn request_limit_is_enforced_while_reading() {
        let mut at_limit = RESOLVE.as_bytes().to_vec();
        at_limit.resize(MAX_REQUEST_BYTES, b' ');
        let parsed: ResolveRequest = decode_request(Cursor::new(&at_limit)).unwrap();
        assert_eq!(parsed.reference, "pdb:1ABC");
        let mut over_limit = at_limit;
        over_limit.resize(MAX_REQUEST_BYTES * 4, b' ');
        let mut reader = Cursor::new(&over_limit);
        let error = decode_request::<ResolveRequest>(&mut reader).unwrap_err();
        assert!(error.contains("64 KiB"));
        assert_eq!(reader.position(), (MAX_REQUEST_BYTES + 1) as u64);
    }

    #[test]
    fn typed_json_rejects_malformed_missing_unknown_duplicate_and_trailing_input() {
        for invalid in [
            "",
            "{",
            "null",
            "[]",
            "{}",
            r#"{"reference":"pdb:1ABC"}"#,
            r#"{"reference":42,"representation":"mmcif"}"#,
            r#"{"reference":"pdb:1ABC","representation":"mmcif","extra":true}"#,
            r#"{"reference":"pdb:1ABC","reference":"pdb:2ABC","representation":"mmcif"}"#,
            r#"{"reference":"pdb:1ABC","representation":"mmcif"} {}"#,
        ] {
            assert!(
                decode_request::<ResolveRequest>(invalid.as_bytes()).is_err(),
                "accepted {invalid}"
            );
        }
        let invalid_declaration = br#"{"source_path":"pkg","requested_ref":"pdb:1ABC","canonical_ref":"pdb:1ABC","declaration":{"provider":"pdb","scope":"entry","representations":{"mmcif":["entry.cif"]},"invented":true}}"#;
        assert!(decode_request::<RegisterRequest>(&invalid_declaration[..]).is_err());
    }
}
