//! Transport-independent, bounded local datasets executed by Rust Polars.
//! Dataset handles are session capabilities, never biological accessions.
use base64::{engine::general_purpose::STANDARD, Engine};
use biov_identifiers::IdentifierRef;
use polars::prelude::*;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    io::{Cursor, Read, Seek, SeekFrom, Write},
    path::{Component, Path, PathBuf},
};
use uuid::Uuid;

pub const MAX_INPUT_BYTES: u64 = 16 * 1024 * 1024;
pub const MAX_PREVIEW_ROWS: usize = 50;
const MAX_DATASETS: usize = 16;
const MAX_MEMORY_BYTES: usize = 64 * 1024 * 1024;
const MAX_COLUMNS: usize = 64;
const PREVIEW_BUDGET: usize = 24 * 1024;
const MAX_CHUNK_BYTES: usize = 48 * 1024;

#[derive(Debug, thiserror::Error)]
#[error("{0}")]
pub struct DatasetError(String);
pub type Result<T> = std::result::Result<T, DatasetError>;
fn err(message: impl Into<String>) -> DatasetError {
    DatasetError(message.into().chars().take(512).collect())
}
fn checked<T, E: std::fmt::Display>(result: std::result::Result<T, E>) -> Result<T> {
    result.map_err(|e| err(e.to_string()))
}
fn default_preview() -> usize {
    5
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema, Default)]
#[serde(deny_unknown_fields)]
pub struct ScientificMetadata {
    /// Caller-declared biological accession. Locally validated, not provider-resolved.
    pub identifier: Option<String>,
    /// Caller-declared species; no species is inferred from column names or identifiers.
    pub species: Option<String>,
    pub reference: Option<String>,
    pub coordinates: Option<String>,
    pub units: Option<String>,
}
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ColumnType {
    String,
    Int64,
    Float64,
    Boolean,
}
impl ColumnType {
    fn dtype(&self) -> DataType {
        match self {
            Self::String => DataType::String,
            Self::Int64 => DataType::Int64,
            Self::Float64 => DataType::Float64,
            Self::Boolean => DataType::Boolean,
        }
    }
}
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct OpenRequest {
    pub path: String,
    /// Unspecified columns remain strings to preserve identifier spelling and numeric precision.
    #[serde(default)]
    pub schema: BTreeMap<String, ColumnType>,
    #[serde(default)]
    pub metadata: ScientificMetadata,
    #[serde(default = "default_preview")]
    pub preview_rows: usize,
}
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PreviewRequest {
    pub dataset_id: String,
    #[serde(default = "default_preview")]
    pub preview_rows: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum FilterOp {
    Eq,
    Gt,
    Ge,
    Lt,
    Le,
    IsNull,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Filter {
    pub column: String,
    pub op: FilterOp,
    pub value: Option<Value>,
}
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Sort {
    pub column: String,
    #[serde(default)]
    pub descending: bool,
}
#[derive(Debug, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct QueryRequest {
    pub dataset_id: String,
    pub filter: Option<Filter>,
    pub select: Option<Vec<String>>,
    pub sort: Option<Sort>,
    #[serde(default = "default_preview")]
    pub preview_rows: usize,
}
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ExportRequest {
    pub dataset_id: String,
}
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReadArtifactRequest {
    pub artifact_id: String,
    #[serde(default)]
    pub offset: u64,
    pub max_bytes: usize,
}
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReleaseRequest {
    pub dataset_id: String,
}

struct Dataset {
    frame: DataFrame,
    provenance: Value,
    metadata: ScientificMetadata,
}
struct Artifact {
    path: PathBuf,
    sha256: String,
    bytes: u64,
}
pub struct DatasetStore {
    data_root: PathBuf,
    output_root: PathBuf,
    datasets: HashMap<String, Dataset>,
    artifacts: HashMap<String, Artifact>,
}

impl DatasetStore {
    /// Roots must already exist and be trusted, operator-owned directories.
    pub fn new(data_root: &Path, output_root: &Path) -> Result<Self> {
        let data_root = checked(data_root.canonicalize())?;
        let output_root = checked(output_root.canonicalize())?;
        if data_root.to_str().is_none() || output_root.to_str().is_none() {
            return Err(err("roots must have UTF-8 paths"));
        }
        if !data_root.is_dir() || !output_root.is_dir() {
            return Err(err("roots must be existing directories"));
        }
        Ok(Self {
            data_root,
            output_root,
            datasets: HashMap::new(),
            artifacts: HashMap::new(),
        })
    }
    fn get(&self, id: &str) -> Result<&Dataset> {
        self.datasets
            .get(id)
            .ok_or_else(|| err("unknown dataset_id (handles last for this server session only)"))
    }
    fn insert(
        &mut self,
        frame: DataFrame,
        provenance: Value,
        metadata: ScientificMetadata,
        rows: usize,
    ) -> Result<Value> {
        validate_frame(&frame)?;
        if checked(serde_json::to_vec(&provenance))?.len() > 8 * 1024 {
            return Err(err("provenance exceeds bounded 8 KiB record limit"));
        }
        if self.datasets.len() >= MAX_DATASETS {
            return Err(err("dataset limit reached; release a dataset first"));
        }
        let retained: usize = self
            .datasets
            .values()
            .map(|d| memory_charge(&d.frame))
            .sum();
        if memory_charge(&frame).saturating_add(retained) > MAX_MEMORY_BYTES {
            return Err(err("retained dataset memory limit exceeded"));
        }
        let id = format!("ds_{}", Uuid::new_v4().simple());
        let dataset = Dataset {
            frame,
            provenance,
            metadata,
        };
        let response = summarize(&id, &dataset, rows)?;
        self.datasets.insert(id, dataset);
        Ok(response)
    }
    /// Reads one immutable in-memory byte snapshot before parsing or hashing.
    /// Local roots are a capability boundary, not a sandbox for hostile local processes.
    pub fn open(&mut self, request: OpenRequest) -> Result<Value> {
        check_rows(request.preview_rows)?;
        validate_metadata(&request.metadata)?;
        let relative = Path::new(&request.path);
        if request.path.len() > 1024
            || relative.as_os_str().is_empty()
            || relative
                .components()
                .any(|c| !matches!(c, Component::Normal(_)))
        {
            return Err(err("path must be a relative file path without traversal"));
        }
        let path = checked(self.data_root.join(relative).canonicalize())?;
        if !path.starts_with(&self.data_root) {
            return Err(err("path is outside configured data root"));
        }
        if !path.is_file() {
            return Err(err("input must be a regular file"));
        }
        if path.extension().and_then(|s| s.to_str()) != Some("csv") {
            return Err(err("this slice opens UTF-8 CSV files only"));
        }
        let file = checked(File::open(&path))?;
        if checked(file.metadata())?.len() > MAX_INPUT_BYTES {
            return Err(err("input exceeds 16 MiB limit"));
        }
        let mut bytes = Vec::new();
        checked(file.take(MAX_INPUT_BYTES + 1).read_to_end(&mut bytes))?;
        if bytes.len() as u64 > MAX_INPUT_BYTES {
            return Err(err("input exceeds 16 MiB limit"));
        }
        let digest = format!("{:x}", Sha256::digest(&bytes));
        let byte_len = bytes.len();
        checked(std::str::from_utf8(&bytes))?;
        let mut csv = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(bytes.as_slice());
        let headers = checked(csv.headers())?.clone();
        if headers.is_empty() || headers.len() > MAX_COLUMNS {
            return Err(err("table must have 1..64 columns"));
        }
        if headers
            .iter()
            .any(|name| name.is_empty() || name.len() > 128 || name.chars().any(char::is_control))
        {
            return Err(err(
                "column names must be 1..128 bytes without control characters",
            ));
        }
        let unique: std::collections::HashSet<_> = headers.iter().collect();
        if unique.len() != headers.len() {
            return Err(err("duplicate CSV column names are not permitted"));
        }
        if request
            .schema
            .keys()
            .any(|name| !unique.contains(name.as_str()))
        {
            return Err(err("schema names must match existing CSV columns"));
        }
        // Reject ragged rows before Polars rather than silently padding or truncating.
        for (index, row) in csv.records().enumerate() {
            checked(row)?;
            if (index + 1).saturating_mul(headers.len()).saturating_mul(17) > MAX_MEMORY_BYTES {
                return Err(err("input row/column allocation exceeds dataset budget"));
            }
        }
        let mut schema = Schema::with_capacity(headers.len());
        for name in &headers {
            schema.with_column(
                name.into(),
                request
                    .schema
                    .get(name)
                    .map(ColumnType::dtype)
                    .unwrap_or(DataType::String),
            );
        }
        let frame = checked(
            CsvReadOptions::default()
                .with_has_header(true)
                .with_schema(Some(std::sync::Arc::new(schema)))
                .into_reader_with_file_handle(Cursor::new(bytes))
                .finish(),
        )?;
        let provenance = json!({"source": request.path, "source_bytes": byte_len, "source_sha256": digest, "input_consistency": "parsed_same_in_memory_snapshot_as_digest", "operations": [], "csv_schema_policy": "strings_unless_explicitly_typed", "declared_schema": request.schema});
        self.insert(frame, provenance, request.metadata, request.preview_rows)
    }
    pub fn preview(&mut self, request: PreviewRequest) -> Result<Value> {
        summarize(
            &request.dataset_id,
            self.get(&request.dataset_id)?,
            request.preview_rows,
        )
    }
    /// Filter, then stable sort, then projection; no implicit row limit.
    pub fn query(&mut self, request: QueryRequest) -> Result<Value> {
        check_rows(request.preview_rows)?;
        let source = self.get(&request.dataset_id)?;
        let mut frame = source.frame.clone();
        if let Some(filter) = &request.filter {
            frame = checked(frame.filter(&filter_mask(&frame, filter)?))?;
        }
        if let Some(sort) = &request.sort {
            frame = checked(
                frame.sort(
                    [sort.column.as_str()],
                    SortMultipleOptions::default()
                        .with_order_descending(sort.descending)
                        .with_maintain_order(true)
                        .with_nulls_last(true),
                ),
            )?;
        }
        if let Some(columns) = &request.select {
            if columns.is_empty() || columns.len() > MAX_COLUMNS {
                return Err(err("select must contain 1..64 unique existing columns"));
            }
            let unique: std::collections::HashSet<_> = columns.iter().collect();
            if unique.len() != columns.len() {
                return Err(err("select columns must be unique"));
            }
            frame = checked(frame.select(columns.iter().map(String::as_str)))?;
        }
        let mut provenance = source.provenance.clone();
        let operation = json!({"parent_dataset_id": request.dataset_id, "filter": request.filter, "sort": request.sort, "select": request.select, "execution_order": ["filter", "stable_sort_nulls_last", "select"]});
        let operations = provenance["operations"]
            .as_array_mut()
            .ok_or_else(|| err("invalid provenance"))?;
        if operations.len() >= 32 {
            return Err(err(
                "maximum 32 derivation steps reached; export and start a new analysis",
            ));
        }
        operations.push(operation);
        if checked(serde_json::to_vec(&provenance))?.len() > 8 * 1024 {
            return Err(err("provenance exceeds bounded 8 KiB record limit"));
        }
        let metadata = source.metadata.clone();
        self.insert(frame, provenance, metadata, request.preview_rows)
    }
    /// Saves complete data and its record; unique names and create-new persistence prevent overwrites.
    pub fn export(&mut self, request: ExportRequest) -> Result<Value> {
        if self.artifacts.len() >= 64 {
            return Err(err("session artifact limit reached"));
        }
        let dataset = self.get(&request.dataset_id)?;
        let artifact_id = format!("artifact_{}", Uuid::new_v4().simple());
        let filename = format!("{artifact_id}.arrow");
        let path = self.output_root.join(&filename);
        let mut temporary = checked(tempfile::NamedTempFile::new_in(&self.output_root))?;
        let mut frame = dataset.frame.clone();
        checked(IpcWriter::new(temporary.as_file_mut()).finish(&mut frame))?;
        checked(temporary.as_file_mut().sync_all())?;
        let bytes = checked(temporary.as_file().metadata())?.len();
        checked(temporary.as_file_mut().seek(SeekFrom::Start(0)))?;
        let mut hasher = Sha256::new();
        checked(std::io::copy(temporary.as_file_mut(), &mut hasher))?;
        let digest = format!("{:x}", hasher.finalize());
        let record = json!({"record_version": 1, "artifact_id": artifact_id, "format": "arrow_ipc", "file": filename, "bytes": bytes, "sha256": digest, "row_count": frame.height(), "schema": schema(&frame), "scientific_metadata": dataset.metadata, "metadata_status": "caller_declared_or_unknown_not_provider_verified", "provenance": dataset.provenance, "software": {"biov": env!("CARGO_PKG_VERSION"), "polars": "0.51.0"}});
        let record_path = self.output_root.join(format!("{artifact_id}.json"));
        let mut record_file = checked(tempfile::NamedTempFile::new_in(&self.output_root))?;
        checked(record_file.write_all(&checked(serde_json::to_vec_pretty(&record))?))?;
        checked(record_file.as_file_mut().sync_all())?;
        checked(temporary.persist_noclobber(&path))?;
        // If the record write fails, the complete IPC stays independently usable but no success is reported.
        checked(record_file.persist_noclobber(&record_path))?;
        self.artifacts.insert(
            artifact_id.clone(),
            Artifact {
                path: path.clone(),
                sha256: digest.clone(),
                bytes,
            },
        );
        Ok(
            json!({"artifact_id": artifact_id, "format": "arrow_ipc", "bytes": bytes, "sha256": digest, "execution_host_path": path, "record_path": record_path, "record": record, "retrieval": {"tool": "dataset_read_artifact", "max_chunk_bytes": MAX_CHUNK_BYTES, "session_scoped": true}, "retention": "files persist; never automatically deleted"}),
        )
    }
    pub fn read_artifact(&mut self, request: ReadArtifactRequest) -> Result<Value> {
        if request.max_bytes == 0 || request.max_bytes > MAX_CHUNK_BYTES {
            return Err(err("max_bytes must be between 1 and 49152"));
        }
        let artifact = self
            .artifacts
            .get(&request.artifact_id)
            .ok_or_else(|| err("unknown artifact_id"))?;
        if request.offset > artifact.bytes {
            return Err(err("offset exceeds artifact size"));
        }
        let file = checked(File::open(&artifact.path))?;
        if checked(file.metadata())?.len() != artifact.bytes {
            return Err(err("artifact changed since export"));
        }
        let mut bytes = Vec::new();
        checked(file.take(artifact.bytes + 1).read_to_end(&mut bytes))?;
        if bytes.len() as u64 != artifact.bytes
            || format!("{:x}", Sha256::digest(&bytes)) != artifact.sha256
        {
            return Err(err(
                "artifact changed since export; refusing mismatched content",
            ));
        }
        let start = request.offset as usize;
        let end = start.saturating_add(request.max_bytes).min(bytes.len());
        Ok(
            json!({"artifact_id": request.artifact_id, "offset": start, "next_offset": end, "total_bytes": bytes.len(), "sha256": artifact.sha256, "encoding": "base64", "data": STANDARD.encode(&bytes[start..end]), "eof": end == bytes.len()}),
        )
    }
    pub fn release(&mut self, request: ReleaseRequest) -> Result<Value> {
        if self.datasets.remove(&request.dataset_id).is_none() {
            return Err(err("unknown dataset_id"));
        }
        Ok(json!({"released": request.dataset_id, "exported_files_preserved": true}))
    }
}

// Polars 0.51 string estimates omit the 16-byte view array. Charge a conservative
// 17 bytes per cell in addition to payload estimates (including null validity).
fn memory_charge(frame: &DataFrame) -> usize {
    frame.estimated_size().saturating_add(
        frame
            .height()
            .saturating_mul(frame.width())
            .saturating_mul(17),
    )
}
fn validate_metadata(metadata: &ScientificMetadata) -> Result<()> {
    for value in [
        &metadata.identifier,
        &metadata.species,
        &metadata.reference,
        &metadata.coordinates,
        &metadata.units,
    ]
    .into_iter()
    .flatten()
    {
        if value.len() > 256 || value.chars().any(char::is_control) {
            return Err(err(
                "metadata values must be at most 256 bytes without control characters",
            ));
        }
    }
    if let Some(identifier) = &metadata.identifier {
        checked(IdentifierRef::parse(identifier))?;
    }
    Ok(())
}
fn validate_frame(frame: &DataFrame) -> Result<()> {
    if frame.width() == 0 || frame.width() > MAX_COLUMNS {
        return Err(err("table must have 1..64 columns"));
    }
    for column in frame.get_columns() {
        if column.dtype() == &DataType::Float64
            && checked(column.f64())?
                .into_iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(err("non-finite float64 values are not supported; preserve them as strings or clean explicitly"));
        }
        if column.name().is_empty()
            || column.name().len() > 128
            || column.name().chars().any(char::is_control)
        {
            return Err(err(
                "column names must be 1..128 bytes without control characters",
            ));
        }
        if !matches!(
            column.dtype(),
            DataType::String | DataType::Int64 | DataType::Float64 | DataType::Boolean
        ) {
            return Err(err(
                "unsupported schema: expected string, int64, float64, or boolean columns",
            ));
        }
    }
    Ok(())
}
fn check_rows(rows: usize) -> Result<()> {
    if rows > MAX_PREVIEW_ROWS {
        Err(err("preview_rows must be at most 50"))
    } else {
        Ok(())
    }
}
fn schema(frame: &DataFrame) -> Value {
    Value::Array(frame.get_columns().iter().map(|column| json!({"name": column.name().as_str(), "dtype": column.dtype().to_string()})).collect())
}
fn summarize(id: &str, dataset: &Dataset, count: usize) -> Result<Value> {
    check_rows(count)?;
    let frame = &dataset.frame;
    let mut rows = Vec::new();
    let mut remaining = PREVIEW_BUDGET;
    let mut truncated_cells = 0;
    for index in 0..count.min(frame.height()) {
        let mut row = Vec::new();
        let mut row_truncated = 0;
        for column in frame.get_columns() {
            let value = checked(column.get(index))?;
            row.push(match value {
                AnyValue::Null => Value::Null,
                AnyValue::Int64(v) => json!(v),
                AnyValue::Float64(v) if v.is_finite() => json!(v),
                AnyValue::Boolean(v) => json!(v),
                AnyValue::StringOwned(v) => {
                    let (text, truncated) = bounded_text(v.as_str());
                    row_truncated += usize::from(truncated);
                    json!(text)
                }
                AnyValue::String(v) => {
                    let (text, truncated) = bounded_text(v);
                    row_truncated += usize::from(truncated);
                    json!(text)
                }
                _ => {
                    row_truncated += 1;
                    json!({"non_finite": value.to_string()})
                }
            });
        }
        let size = checked(serde_json::to_vec(&row))?.len();
        if size > remaining {
            break;
        }
        remaining -= size;
        truncated_cells += row_truncated;
        rows.push(row);
    }
    Ok(
        json!({"dataset_id": id, "row_count": frame.height(), "column_count": frame.width(), "schema": schema(frame), "preview": {"rows": rows, "requested_rows": count, "returned_rows": rows.len(), "omitted_rows": frame.height()-rows.len(), "truncated_cells": truncated_cells, "max_cell_bytes": 256, "row_order": "source order unless explicit stable sort", "not_complete_data": rows.len() != frame.height() || truncated_cells > 0}, "scientific_metadata": dataset.metadata, "metadata_status": "caller_declared_or_unknown_not_provider_verified", "source_sha256": dataset.provenance["source_sha256"], "handle_lifetime": "server_session", "last_operation": dataset.provenance["operations"].as_array().and_then(|ops| ops.last()), "engine": "rust_polars"}),
    )
}
fn bounded_text(value: &str) -> (String, bool) {
    if value.len() <= 256 {
        return (value.into(), false);
    }
    let mut end = 253;
    while !value.is_char_boundary(end) {
        end -= 1;
    }
    (format!("{}...", &value[..end]), true)
}
fn filter_mask(frame: &DataFrame, filter: &Filter) -> Result<BooleanChunked> {
    let column = checked(frame.column(&filter.column))?;
    if matches!(filter.op, FilterOp::IsNull) {
        if filter.value.is_some() {
            return Err(err("is_null must omit value"));
        }
        return Ok(column.is_null());
    }
    let value = filter
        .value
        .as_ref()
        .ok_or_else(|| err("filter requires a typed non-null value"))?;
    macro_rules! comparison {
        ($column:expr, $value:expr) => {
            match filter.op {
                FilterOp::Eq => $column.equal($value),
                FilterOp::Gt => $column.gt($value),
                FilterOp::Ge => $column.gt_eq($value),
                FilterOp::Lt => $column.lt($value),
                FilterOp::Le => $column.lt_eq($value),
                FilterOp::IsNull => unreachable!(),
            }
        };
    }
    Ok(match column.dtype() {
        DataType::Int64 => comparison!(
            checked(column.i64())?,
            value
                .as_i64()
                .ok_or_else(|| err("int64 filter requires an exact int64 value"))?
        ),
        DataType::Float64 => comparison!(
            checked(column.f64())?,
            value
                .as_f64()
                .filter(|v| v.is_finite())
                .ok_or_else(|| err("float64 filter requires a finite numeric value"))?
        ),
        DataType::String => comparison!(
            checked(column.str())?,
            value
                .as_str()
                .ok_or_else(|| err("string filter requires a string value"))?
        ),
        DataType::Boolean => {
            if !matches!(filter.op, FilterOp::Eq) {
                return Err(err("boolean filters support eq or is_null only"));
            }
            let target = value
                .as_bool()
                .ok_or_else(|| err("boolean filter requires a boolean value"))?;
            let boolean = checked(column.bool())?;
            if target {
                boolean.clone()
            } else {
                !boolean
            }
        }
        _ => return Err(err("unsupported filter dtype")),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn setup(csv: &str) -> (tempfile::TempDir, DatasetStore) {
        let dir = tempfile::tempdir().unwrap();
        fs::create_dir(dir.path().join("inputs")).unwrap();
        fs::create_dir(dir.path().join("outputs")).unwrap();
        fs::write(dir.path().join("inputs/table.csv"), csv).unwrap();
        let store =
            DatasetStore::new(&dir.path().join("inputs"), &dir.path().join("outputs")).unwrap();
        (dir, store)
    }
    fn open(store: &mut DatasetStore, preview_rows: usize) -> Value {
        store
            .open(OpenRequest {
                schema: if std::fs::read_to_string(store.data_root.join("table.csv"))
                    .unwrap()
                    .lines()
                    .next()
                    .unwrap()
                    .split(',')
                    .any(|name| name == "count")
                {
                    BTreeMap::from([("count".into(), ColumnType::Int64)])
                } else {
                    BTreeMap::new()
                },
                path: "table.csv".into(),
                metadata: ScientificMetadata::default(),
                preview_rows,
            })
            .unwrap()
    }
    #[test]
    fn default_strings_preserve_identifier_spelling_and_precision() {
        let (_dir, mut store) = setup("id,value\n00123,9007199254740993\n00456,1.5\n");
        let result = open(&mut store, 5);
        assert_eq!(
            result["preview"]["rows"],
            json!([["00123", "9007199254740993"], ["00456", "1.5"]])
        );
    }
    #[test]
    fn invalid_typed_schema_and_nonfinite_numbers_are_rejected() {
        for csv in ["value\nNaN\n", "value\ninf\n", "value\nabc\n"] {
            let (_dir, mut store) = setup(csv);
            assert!(store
                .open(OpenRequest {
                    path: "table.csv".into(),
                    schema: BTreeMap::from([("value".into(), ColumnType::Float64)]),
                    metadata: ScientificMetadata::default(),
                    preview_rows: 1
                })
                .is_err());
        }
        for csv in ["a,a\n1,2\n", "a,b\n1,2,3\n"] {
            let (_dir, mut store) = setup(csv);
            assert!(store
                .open(OpenRequest {
                    path: "table.csv".into(),
                    schema: BTreeMap::new(),
                    metadata: ScientificMetadata::default(),
                    preview_rows: 1
                })
                .is_err());
        }
    }
    #[test]
    fn complete_data_is_filtered_and_ipc_roundtrips() {
        let (_dir, mut store) = setup("sample,count\na,1\nb,30\nc,20\nd,30\n");
        let initial = open(&mut store, 1);
        assert_eq!(initial["row_count"], 4);
        assert_eq!(initial["scientific_metadata"]["species"], Value::Null);
        let result = store.query(serde_json::from_value(json!({"dataset_id": initial["dataset_id"], "filter": {"column":"count", "op":"ge", "value":20}, "sort":{"column":"count","descending":true}, "select":["sample","count"], "preview_rows":1})).unwrap()).unwrap();
        assert_eq!(result["row_count"], 3);
        assert_eq!(result["preview"]["rows"], json!([["b", 30]]));
        let exported = store
            .export(ExportRequest {
                dataset_id: result["dataset_id"].as_str().unwrap().into(),
            })
            .unwrap();
        let frame =
            IpcReader::new(File::open(exported["execution_host_path"].as_str().unwrap()).unwrap())
                .finish()
                .unwrap();
        assert_eq!(frame.height(), 3);
        assert_eq!(
            frame
                .column("sample")
                .unwrap()
                .str()
                .unwrap()
                .into_iter()
                .collect::<Vec<_>>(),
            vec![Some("b"), Some("d"), Some("c")]
        );
        assert_eq!(
            exported["record"]["provenance"]["operations"]
                .as_array()
                .unwrap()
                .len(),
            1
        );
        let mut all = Vec::new();
        loop {
            let chunk = store
                .read_artifact(ReadArtifactRequest {
                    artifact_id: exported["artifact_id"].as_str().unwrap().into(),
                    offset: all.len() as u64,
                    max_bytes: 17,
                })
                .unwrap();
            all.extend(STANDARD.decode(chunk["data"].as_str().unwrap()).unwrap());
            if chunk["eof"] == true {
                break;
            }
        }
        assert_eq!(format!("{:x}", Sha256::digest(&all)), exported["sha256"]);
        assert_eq!(
            IpcReader::new(Cursor::new(all)).finish().unwrap().height(),
            3
        );
        assert_eq!(
            initial["source_sha256"],
            exported["record"]["provenance"]["source_sha256"]
        );
    }
    #[test]
    fn exact_integer_filters_nulls_and_stable_sort() {
        let (_dir, mut store) = setup("label,count\na,9007199254740992\nb,9007199254740993\nc,\n");
        let initial = open(&mut store, 5);
        let id = initial["dataset_id"].as_str().unwrap();
        let result = store.query(serde_json::from_value(json!({"dataset_id":id,"filter":{"column":"count","op":"eq","value":9007199254740993_i64}})).unwrap()).unwrap();
        assert_eq!(
            result["preview"]["rows"],
            json!([["b", 9007199254740993_i64]])
        );
        let nulls = store
            .query(
                serde_json::from_value(
                    json!({"dataset_id":id,"filter":{"column":"count","op":"is_null"}}),
                )
                .unwrap(),
            )
            .unwrap();
        assert_eq!(nulls["preview"]["rows"], json!([["c", null]]));
        assert!(store
            .query(
                serde_json::from_value(
                    json!({"dataset_id":id,"filter":{"column":"count","op":"eq","value":1.5}})
                )
                .unwrap()
            )
            .is_err());
    }
    #[test]
    fn path_limits_and_metadata_are_explicit() {
        let (dir, mut store) = setup("x\n1\n");
        for path in ["../table.csv", "/etc/passwd", "missing.csv", ""] {
            assert!(store
                .open(OpenRequest {
                    schema: BTreeMap::new(),
                    path: path.into(),
                    metadata: ScientificMetadata::default(),
                    preview_rows: 1
                })
                .is_err());
        }
        assert!(store
            .open(OpenRequest {
                schema: BTreeMap::new(),
                path: "table.csv".into(),
                metadata: ScientificMetadata {
                    identifier: Some("ds_fake".into()),
                    ..Default::default()
                },
                preview_rows: 1
            })
            .is_err());
        assert!(store
            .open(OpenRequest {
                schema: BTreeMap::new(),
                path: "table.csv".into(),
                metadata: ScientificMetadata::default(),
                preview_rows: 51
            })
            .is_err());
        #[cfg(unix)]
        {
            fs::write(dir.path().join("outside.csv"), "x\n2\n").unwrap();
            std::os::unix::fs::symlink(
                dir.path().join("outside.csv"),
                dir.path().join("inputs/escape.csv"),
            )
            .unwrap();
            assert!(store
                .open(OpenRequest {
                    schema: BTreeMap::new(),
                    path: "escape.csv".into(),
                    metadata: ScientificMetadata::default(),
                    preview_rows: 1
                })
                .is_err());
        }
        let file = File::create(dir.path().join("inputs/large.csv")).unwrap();
        file.set_len(MAX_INPUT_BYTES + 1).unwrap();
        assert!(store
            .open(OpenRequest {
                schema: BTreeMap::new(),
                path: "large.csv".into(),
                metadata: ScientificMetadata::default(),
                preview_rows: 1
            })
            .is_err());
    }
    #[test]
    fn release_invalid_handles_and_artifact_tampering() {
        let (_dir, mut store) = setup("x\n1\n");
        let value = open(&mut store, 1);
        let id = value["dataset_id"].as_str().unwrap().to_owned();
        let first = store
            .export(ExportRequest {
                dataset_id: id.clone(),
            })
            .unwrap();
        let second = store
            .export(ExportRequest {
                dataset_id: id.clone(),
            })
            .unwrap();
        assert_ne!(first["artifact_id"], second["artifact_id"]);
        store
            .release(ReleaseRequest {
                dataset_id: id.clone(),
            })
            .unwrap();
        assert!(store
            .preview(PreviewRequest {
                dataset_id: id,
                preview_rows: 1
            })
            .is_err());
        assert!(Path::new(first["execution_host_path"].as_str().unwrap()).exists());
        let artifact_id = first["artifact_id"].as_str().unwrap().to_owned();
        for max_bytes in [0, MAX_CHUNK_BYTES + 1] {
            assert!(store
                .read_artifact(ReadArtifactRequest {
                    artifact_id: artifact_id.clone(),
                    offset: 0,
                    max_bytes
                })
                .is_err());
        }
        fs::write(first["execution_host_path"].as_str().unwrap(), "changed").unwrap();
        assert!(store
            .read_artifact(ReadArtifactRequest {
                artifact_id,
                offset: 0,
                max_bytes: 10
            })
            .is_err());
    }
    #[test]
    fn unicode_preview_is_bounded_without_changing_full_data() {
        let text = "基因".repeat(10000);
        let (_dir, mut store) = setup(&format!("description\n{text}\n"));
        let preview = open(&mut store, 1);
        assert_eq!(preview["preview"]["truncated_cells"], 1);
        assert!(preview["preview"]["rows"][0][0].as_str().unwrap().len() <= 256);
        let result = store
            .export(ExportRequest {
                dataset_id: preview["dataset_id"].as_str().unwrap().into(),
            })
            .unwrap();
        let frame =
            IpcReader::new(File::open(result["execution_host_path"].as_str().unwrap()).unwrap())
                .finish()
                .unwrap();
        assert_eq!(
            frame.column("description").unwrap().str().unwrap().get(0),
            Some(text.as_str())
        );
    }
    #[test]
    fn schema_and_session_count_are_bounded() {
        let headers = (0..65)
            .map(|i| format!("col{i}"))
            .collect::<Vec<_>>()
            .join(",");
        let (_dir, mut wide) = setup(&format!("{headers}\n{}\n", vec!["1"; 65].join(",")));
        assert!(wide
            .open(OpenRequest {
                schema: BTreeMap::new(),
                path: "table.csv".into(),
                metadata: ScientificMetadata::default(),
                preview_rows: 1
            })
            .is_err());
        let (_dir, mut store) = setup("x\n1\n");
        for _ in 0..MAX_DATASETS {
            open(&mut store, 0);
        }
        assert!(store
            .open(OpenRequest {
                schema: BTreeMap::new(),
                path: "table.csv".into(),
                metadata: ScientificMetadata::default(),
                preview_rows: 0
            })
            .is_err());
    }
}
