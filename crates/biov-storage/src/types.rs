use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Resource limits apply to metadata, traversal and response cardinality, not
/// aggregate native file bytes. All file contents are copied/hashed in chunks.
pub const MAX_ENTRIES: usize = 100_000;
pub const MAX_METADATA_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_PATH_BYTES: usize = 4096;
pub const MAX_REPRESENTATIONS: usize = 128;
pub const MAX_RESOLVED_FILES: usize = 256;
pub const MAX_SCAN_CANDIDATES: usize = 1024;
pub const MAX_AMBIGUOUS_CANDIDATES: usize = 64;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "provider", rename_all = "snake_case", deny_unknown_fields)]
pub enum NativeDeclaration {
    /// An extracted, materialized NCBI Datasets RefSeq GCF package.
    Refseq,
    /// Caller declarations only; no biological identity or format validation.
    Pdb {
        /// `entry` or `assembly:<positive decimal ID>`.
        scope: String,
        /// Representation names to nonempty vectors of source-relative paths.
        #[serde(deserialize_with = "unique_representations")]
        representations: BTreeMap<String, Vec<String>>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RegisterRequest {
    pub source_path: String,
    pub requested_ref: String,
    pub canonical_ref: String,
    pub declaration: NativeDeclaration,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ResolveRequest {
    pub reference: String,
    pub representation: String,
    pub snapshot_id: Option<String>,
    /// Defaults to `assembly` for RefSeq and `entry` for PDB.
    pub scope: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum EntryKind {
    Directory,
    File,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct InventoryEntry {
    /// UTF-8, slash-separated, source-relative path.
    pub path: String,
    pub kind: EntryKind,
    /// Zero for directories; bytes saved on disk (including gzip if retained).
    pub bytes: u64,
    /// Lowercase hex SHA-256 of file bytes; null for directories.
    pub sha256: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Validation {
    /// Stable validation contract identifier, checked against repeated native validation.
    pub method: String,
    /// Historical, untrusted explanatory prose; not an integrity or authentication claim.
    pub limits: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Receipt {
    pub schema_version: u32,
    pub snapshot_id: String,
    pub source_content_sha256: String,
    /// SHA-256 of the companion README bytes, independent of future rendering.
    pub readme_sha256: String,
    /// First local registration request; subsequent reuse does not rewrite it.
    pub requested_ref: String,
    pub canonical_ref: String,
    pub scope: String,
    pub declaration: NativeDeclaration,
    /// Native relative paths are authoritative; the wrapper adds navigation.
    pub native_metadata: Vec<String>,
    #[serde(deserialize_with = "unique_representations")]
    pub representations: BTreeMap<String, Vec<String>>,
    pub unavailable_representations: Vec<String>,
    pub inventory: Vec<InventoryEntry>,
    pub validation: Validation,
    pub origin: String,
    pub registered_at_unix_seconds: u64,
    /// Unknown acquisition facts stay null. Registration time is not retrieval.
    pub acquired_at: Option<String>,
    pub source_url: Option<String>,
    pub acquisition_tool: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct SnapshotSummary {
    pub snapshot_id: String,
    pub snapshot_path: String,
    pub receipt_path: String,
    pub source_content_sha256: String,
    pub requested_ref: String,
    pub canonical_ref: String,
    pub scope: String,
    pub file_count: usize,
    pub total_bytes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct Registration {
    pub reused: bool,
    /// This request, distinct from the first request in the immutable receipt.
    pub requested_ref: String,
    pub snapshot: SnapshotSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct ResolvedFile {
    /// Portable path relative to the configured store root.
    pub relative_path: String,
    /// Convenience only; valid on the machine performing resolution.
    pub execution_host_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum Resolution {
    Ready {
        snapshot: SnapshotSummary,
        paths: Vec<ResolvedFile>,
    },
    Miss,
    Ambiguous {
        candidates: Vec<SnapshotSummary>,
    },
    Unavailable {
        snapshot: SnapshotSummary,
        available_representations: Vec<String>,
    },
    Corrupt {
        snapshot_id: String,
        detail: String,
    },
}

#[derive(Debug, thiserror::Error)]
pub enum StorageError {
    #[error("invalid storage request: {0}")]
    InvalidInput(String),
    #[error("native package validation failed: {0}")]
    InvalidPackage(String),
    #[error("storage metadata or output limit exceeded: {0}")]
    Limit(String),
    #[error("snapshot is corrupt: {0}")]
    Corrupt(String),
    #[error("identical source content has a different declared scope or representation mapping")]
    DeclarationConflict,
    #[error("storage filesystem operation failed ({operation}): {kind:?}")]
    Io {
        operation: &'static str,
        kind: std::io::ErrorKind,
    },
    #[error("atomic no-replace directory publication is unsupported on this platform/filesystem")]
    UnsupportedPublication,
}

pub(crate) fn io(operation: &'static str, error: std::io::Error) -> StorageError {
    StorageError::Io {
        operation,
        kind: error.kind(),
    }
}

// A representation map is an ordered-path declaration, not a last-key-wins
// configuration. Reject duplicate JSON keys before they can discard alternatives.
fn unique_representations<'de, D>(
    deserializer: D,
) -> Result<BTreeMap<String, Vec<String>>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    struct UniqueMap;
    impl<'de> serde::de::Visitor<'de> for UniqueMap {
        type Value = BTreeMap<String, Vec<String>>;
        fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            formatter.write_str("a representation object with unique keys and bounded path vectors")
        }
        fn visit_map<A>(self, mut map: A) -> Result<Self::Value, A::Error>
        where
            A: serde::de::MapAccess<'de>,
        {
            let mut entries = BTreeMap::new();
            while let Some((key, paths)) = map.next_entry::<String, Vec<String>>()? {
                if entries.contains_key(&key) {
                    return Err(serde::de::Error::custom("duplicate representation key"));
                }
                if entries.len() >= MAX_REPRESENTATIONS || paths.len() > MAX_RESOLVED_FILES {
                    return Err(serde::de::Error::custom(
                        "representation mapping exceeds count bounds",
                    ));
                }
                entries.insert(key, paths);
            }
            Ok(entries)
        }
    }
    deserializer.deserialize_map(UniqueMap)
}
