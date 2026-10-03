use biov_storage::ResolvedFile;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// Physical line bytes, including LF or CRLF. No whole sequence is retained.
pub const MAX_LINE_BYTES: usize = 1024 * 1024;
pub const MAX_ID_BYTES: usize = 4096;
pub const MAX_RECORDS: usize = 100_000;
/// Applies separately to retained identifier charge and each generated index/table.
pub const MAX_METADATA_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_PROVENANCE_BYTES: usize = 64 * 1024;
pub const IO_BUFFER_BYTES: usize = 64 * 1024;

/// An exact, already registered native RefSeq genome FASTA. No downloads occur.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PrepareFastaRequest {
    /// Exact canonical `refseq.gcf:GCF_000005845.2` form, with biological version.
    pub reference: String,
    /// Exact `sha256-` plus 64 lowercase hexadecimal characters.
    pub snapshot_id: String,
    /// Selected `genome_fasta` file, relative to the native snapshot's `source/`.
    pub source_path: String,
}

/// Compact verified paths. Host paths are conveniences, not portable identities.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
pub struct PreparedFasta {
    pub recipe_id: String,
    pub reused: bool,
    pub reference: String,
    pub snapshot_id: String,
    pub sequence_count: usize,
    pub total_bases: u64,
    pub fasta: ResolvedFile,
    pub fai: ResolvedFile,
    pub dictionary: ResolvedFile,
    pub provenance: ResolvedFile,
    pub readme: ResolvedFile,
}

/// Deterministic recipe inputs, deliberately separate from output checksums.
/// The contract revision MUST change for any output-affecting implementation change.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FastaRecipe {
    pub schema_version: u32,
    pub algorithm: String,
    pub algorithm_revision: u32,
    pub implementation: String,
    pub implementation_version: String,
    pub library: String,
    pub library_version: String,
    pub reference: String,
    pub snapshot_id: String,
    pub source_content_sha256: String,
    pub source_path: String,
    pub source_sha256: String,
    pub source_bytes: u64,
    pub parameters: FastaParameters,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FastaParameters {
    pub input_format: String,
    pub identifier_convention: String,
    pub sequence_byte_policy: String,
    pub dictionary_schema: String,
    pub max_physical_line_bytes: usize,
    pub max_identifier_bytes: usize,
    pub max_records: usize,
    pub max_metadata_bytes: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputIdentity {
    /// Relative to the prepared directory, never an absolute host path.
    pub path: String,
    pub bytes: u64,
    pub sha256: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FastaProvenance {
    pub schema_version: u32,
    pub recipe_id: String,
    pub recipe: FastaRecipe,
    /// Relative to this prepared directory. Copy its exact snapshot dependency too.
    pub source_relative_path: String,
    /// Relative to this prepared directory; the complete native snapshot is required.
    pub source_snapshot_relative_path: String,
    pub sequence_count: usize,
    pub total_bases: u64,
    pub outputs: Vec<OutputIdentity>,
}

#[derive(Debug, thiserror::Error)]
pub enum PreparedError {
    #[error("invalid prepared FASTA request: {0}")]
    InvalidInput(String),
    #[error("unsupported or malformed plain FASTA: {0}")]
    InvalidFasta(String),
    #[error("prepared FASTA limit exceeded: {0}")]
    Limit(String),
    #[error("native input is not ready: {0}")]
    InputUnavailable(String),
    #[error("prepared FASTA or input is corrupt: {0}")]
    Corrupt(String),
    #[error(transparent)]
    Storage(#[from] biov_storage::StorageError),
    #[error("prepared filesystem operation failed ({operation}): {kind:?}")]
    Io {
        operation: &'static str,
        kind: std::io::ErrorKind,
    },
    #[error("atomic no-replace publication is unsupported on this platform/filesystem")]
    UnsupportedPublication,
}

pub(crate) fn io(operation: &'static str, error: std::io::Error) -> PreparedError {
    PreparedError::Io {
        operation,
        kind: error.kind(),
    }
}
pub(crate) fn invalid(message: &str) -> PreparedError {
    PreparedError::InvalidInput(message.into())
}
pub(crate) fn corrupt(message: &str) -> PreparedError {
    PreparedError::Corrupt(message.into())
}
pub(crate) fn limit(message: &str) -> PreparedError {
    PreparedError::Limit(message.into())
}
