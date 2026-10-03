//! Reopen only bounded native BioV exports, not arbitrary Arrow IPC imports.
//! A matching record proves consistency with that record, never its authenticity.
use super::*;
use polars_arrow_format::ipc::{self, planus::ReadAsRoot};
use std::collections::HashSet;

pub(super) const MAX_RECORD_BYTES: u64 = 64 * 1024;
pub(super) const MAX_IPC_BYTES: u64 = 65 * 1024 * 1024;
const MAX_IPC_METADATA: usize = 1024 * 1024;
const MAX_BATCHES: usize = 4096;
const MAX_BUFFERS: usize = 4096;
const METADATA_STATUS: &str = "caller_declared_or_unknown_not_provider_verified";
const CONSISTENCY: &str = "parsed_same_in_memory_snapshot_as_digest";
const CLAIMS: &str = "recorded_claims_not_independently_verified";
const CHECKS: [&str; 4] = ["artifact_bytes", "artifact_sha256", "schema", "row_count"];

#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReopenRequest {
    /// Existing BioV JSON record relative to the configured data root.
    pub record_path: String,
    #[serde(default = "default_preview")]
    pub preview_rows: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct ReopenVerification {
    record_path: String,
    record_sha256: String,
    artifact_sha256: String,
    artifact_bytes: u64,
    record_version: u32,
    checks: Vec<String>,
    input_consistency: String,
    original_provenance: String,
    authenticity: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ExportRecord {
    record_version: u32,
    artifact_id: String,
    format: String,
    file: String,
    bytes: u64,
    sha256: String,
    row_count: usize,
    schema: Vec<RecordedColumn>,
    scientific_metadata: ScientificMetadata,
    metadata_status: String,
    provenance: Provenance,
    software: Software,
    #[serde(default)]
    reopen_verification: Option<ReopenVerification>,
}
#[derive(Deserialize, Serialize, Debug, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
struct RecordedColumn {
    name: String,
    dtype: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Software {
    biov: String,
    polars: String,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Provenance {
    source: String,
    source_bytes: u64,
    source_sha256: String,
    input_consistency: String,
    operations: Vec<Operation>,
    csv_schema_policy: String,
    declared_schema: BTreeMap<String, ColumnType>,
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Operation {
    parent_dataset_id: String,
    filter: Option<Filter>,
    sort: Option<Sort>,
    select: Option<Vec<String>>,
    execution_order: Vec<String>,
}

impl DatasetStore {
    /// Reads the record and IPC as separately bounded snapshots. The IPC bytes
    /// hashed, preflighted, and parsed are the same immutable in-memory bytes.
    pub fn reopen(&mut self, request: ReopenRequest) -> Result<Value> {
        check_rows(request.preview_rows)?;
        if self.datasets.len() >= MAX_DATASETS {
            return Err(err("dataset limit reached; release a dataset first"));
        }
        let relative = relative_path(&request.record_path)?;
        if relative.extension().and_then(|s| s.to_str()) != Some("json") {
            return Err(err("reopen requires a BioV JSON export record"));
        }
        let record_path = confined_file(&self.data_root, relative)?;
        let record_bytes = snapshot(&record_path, MAX_RECORD_BYTES, "export record")?;
        // Deserialize directly into strict structs, including nested provenance,
        // so unknown and duplicate fields are rejected rather than discarded.
        let record: ExportRecord = checked(serde_json::from_slice(&record_bytes))?;
        validate_record(&record)?;
        let filename = relative_path(&record.file)?;
        if filename.components().count() != 1 {
            return Err(err("record IPC file must be a same-directory basename"));
        }
        let record_parent = record_path
            .parent()
            .ok_or_else(|| err("record has no parent directory"))?;
        let ipc_path = confined_file(&self.data_root, &record_parent.join(filename))?;
        let ipc_bytes = snapshot(&ipc_path, MAX_IPC_BYTES, "IPC artifact")?;
        if ipc_bytes.len() as u64 != record.bytes {
            return Err(err("IPC artifact byte size does not match export record"));
        }
        let digest = format!("{:x}", Sha256::digest(&ipc_bytes));
        if digest != record.sha256 {
            return Err(err("IPC artifact SHA-256 does not match export record"));
        }
        let retained = self.datasets.values().fold(0usize, |total, dataset| {
            total.saturating_add(memory_charge(&dataset.frame))
        });
        preflight_ipc(&ipc_bytes, &record.schema, record.row_count, retained)?;
        // No mmap or second file read. Disable automatic rechunk allocations;
        // both per-batch buffers and aggregate rows were bounded above.
        let frame = checked(
            IpcReader::new(Cursor::new(ipc_bytes))
                .set_rechunk(false)
                .finish(),
        )?;
        validate_frame(&frame)?;
        if frame.height() != record.row_count
            || schema(&frame) != checked(serde_json::to_value(&record.schema))?
        {
            return Err(err(
                "parsed IPC schema or row count does not match export record",
            ));
        }
        let verification = ReopenVerification {
            record_path: request.record_path,
            record_sha256: format!("{:x}", Sha256::digest(&record_bytes)),
            artifact_sha256: digest,
            artifact_bytes: record.bytes,
            record_version: record.record_version,
            checks: CHECKS.iter().map(|value| (*value).into()).collect(),
            input_consistency: CONSISTENCY.into(),
            original_provenance: CLAIMS.into(),
            authenticity: "not_established".into(),
        };
        self.insert(
            frame,
            checked(serde_json::to_value(record.provenance))?,
            record.scientific_metadata,
            Some(verification),
            request.preview_rows,
        )
    }
}

fn relative_path(value: &str) -> Result<&Path> {
    let path = Path::new(value);
    if value.is_empty()
        || value.len() > 1024
        || path
            .components()
            .any(|c| !matches!(c, Component::Normal(_)))
    {
        return Err(err("path must be a relative file path without traversal"));
    }
    Ok(path)
}
fn confined_file(root: &Path, relative: &Path) -> Result<PathBuf> {
    let path = checked(root.join(relative).canonicalize())?;
    if !path.starts_with(root) {
        return Err(err("path is outside configured data root"));
    }
    if !path.is_file() {
        return Err(err("input must be a regular file"));
    }
    Ok(path)
}
fn snapshot(path: &Path, cap: u64, label: &str) -> Result<Vec<u8>> {
    let file = checked(File::open(path))?;
    if checked(file.metadata())?.len() > cap {
        return Err(err(format!("{label} exceeds bounded byte limit ({cap})")));
    }
    let mut bytes = Vec::new();
    checked(file.take(cap + 1).read_to_end(&mut bytes))?;
    if bytes.len() as u64 > cap {
        return Err(err(format!("{label} exceeds bounded byte limit ({cap})")));
    }
    Ok(bytes)
}
fn digest_valid(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn session_id(value: &str, prefix: &str) -> bool {
    value.strip_prefix(prefix).is_some_and(|suffix| {
        suffix.len() == 32
            && suffix
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    })
}
fn column_name(value: &str) -> Result<()> {
    if value.is_empty() || value.len() > 128 || value.chars().any(char::is_control) {
        return Err(err(
            "column names must be 1..128 bytes without control characters",
        ));
    }
    Ok(())
}
fn validate_verification(value: &ReopenVerification) -> Result<()> {
    relative_path(&value.record_path)?;
    if !digest_valid(&value.record_sha256)
        || !digest_valid(&value.artifact_sha256)
        || value.artifact_bytes > MAX_IPC_BYTES
        || !matches!(value.record_version, 1 | 2)
        || value.checks != CHECKS
        || value.input_consistency != CONSISTENCY
        || value.original_provenance != CLAIMS
        || value.authenticity != "not_established"
    {
        return Err(err("invalid recorded reopen verification"));
    }
    Ok(())
}
fn validate_record(record: &ExportRecord) -> Result<()> {
    if !matches!(record.record_version, 1 | 2) || record.format != "arrow_ipc" {
        return Err(err("unsupported export record version or format"));
    }
    if record.record_version == 1 && record.reopen_verification.is_some() {
        return Err(err("version 1 records cannot contain reopen verification"));
    }
    if let Some(verification) = &record.reopen_verification {
        validate_verification(verification)?;
    }
    if !session_id(&record.artifact_id, "artifact_")
        || record.file != format!("{}.arrow", record.artifact_id)
    {
        return Err(err(
            "record file and artifact identity must match native export names",
        ));
    }
    if record.bytes == 0 || record.bytes > MAX_IPC_BYTES || !digest_valid(&record.sha256) {
        return Err(err("invalid or oversized IPC content identity"));
    }
    if record.schema.is_empty() || record.schema.len() > MAX_COLUMNS {
        return Err(err("record schema must contain 1..64 columns"));
    }
    let mut names = HashSet::new();
    for column in &record.schema {
        column_name(&column.name)?;
        if !names.insert(&column.name)
            || ![
                DataType::String,
                DataType::Int64,
                DataType::Float64,
                DataType::Boolean,
            ]
            .iter()
            .any(|dtype| dtype.to_string() == column.dtype)
        {
            return Err(err("unsupported or duplicate record schema column"));
        }
    }
    allocation_budget(record.row_count, record.schema.len(), 0)?;
    validate_metadata(&record.scientific_metadata)?;
    if record.metadata_status != METADATA_STATUS
        || record.software.polars != "0.51.0"
        || record.software.biov.is_empty()
        || record.software.biov.len() > 128
        || record.software.biov.chars().any(char::is_control)
    {
        return Err(err("unsupported or malformed record metadata/software"));
    }
    let source = &record.provenance;
    relative_path(&source.source)?;
    if source.source_bytes > MAX_INPUT_BYTES
        || !digest_valid(&source.source_sha256)
        || source.input_consistency != CONSISTENCY
        || source.csv_schema_policy != "strings_unless_explicitly_typed"
        || source.operations.len() > 32
        || source.declared_schema.len() > MAX_COLUMNS
    {
        return Err(err("invalid recorded source provenance"));
    }
    for name in source.declared_schema.keys() {
        column_name(name)?;
    }
    for operation in &source.operations {
        if !session_id(&operation.parent_dataset_id, "ds_")
            || operation.execution_order != ["filter", "stable_sort_nulls_last", "select"]
        {
            return Err(err("invalid recorded operation"));
        }
        if let Some(filter) = &operation.filter {
            column_name(&filter.column)?;
            if matches!(filter.op, FilterOp::IsNull) != filter.value.is_none()
                || filter
                    .value
                    .as_ref()
                    .is_some_and(|v| !(v.is_string() || v.is_boolean() || v.is_number()))
            {
                return Err(err("invalid recorded filter"));
            }
        }
        if let Some(sort) = &operation.sort {
            column_name(&sort.column)?;
        }
        if let Some(select) = &operation.select {
            if select.is_empty()
                || select.len() > MAX_COLUMNS
                || select.iter().collect::<HashSet<_>>().len() != select.len()
            {
                return Err(err("invalid recorded projection"));
            }
            for name in select {
                column_name(name)?;
            }
        }
    }
    if checked(serde_json::to_vec(source))?.len() > 8 * 1024 {
        return Err(err("provenance exceeds bounded 8 KiB record limit"));
    }
    Ok(())
}
fn allocation_budget(rows: usize, columns: usize, retained: usize) -> Result<()> {
    if rows
        .checked_mul(columns)
        .and_then(|cells| cells.checked_mul(17))
        .and_then(|charge| charge.checked_add(retained))
        .is_none_or(|charge| charge > MAX_MEMORY_BYTES)
    {
        return Err(err(
            "IPC row/column allocation exceeds retained dataset budget",
        ));
    }
    Ok(())
}
fn nonnegative(value: i64) -> Result<usize> {
    usize::try_from(value).map_err(|_| err("negative or oversized IPC metadata integer"))
}
fn checked_end(start: usize, length: usize, end: usize) -> Result<usize> {
    start
        .checked_add(length)
        .filter(|value| *value <= end)
        .ok_or_else(|| err("IPC metadata range is outside bounded snapshot"))
}
fn slice(bytes: &[u8], start: usize, length: usize) -> Result<&[u8]> {
    Ok(&bytes[start..checked_end(start, length, bytes.len())?])
}
#[derive(Clone, Copy, PartialEq, Eq)]
enum NativeType {
    StringView,
    LargeString,
    Int64,
    Float64,
    Boolean,
}
impl NativeType {
    fn dtype(self) -> DataType {
        match self {
            Self::StringView | Self::LargeString => DataType::String,
            Self::Int64 => DataType::Int64,
            Self::Float64 => DataType::Float64,
            Self::Boolean => DataType::Boolean,
        }
    }
    fn is_string(self) -> bool {
        matches!(self, Self::StringView | Self::LargeString)
    }
}
fn native_schema(schema: ipc::SchemaRef<'_>) -> Result<Vec<(String, NativeType)>> {
    if checked(schema.endianness())? != ipc::Endianness::Little
        || !cfg!(target_endian = "little")
        || checked(schema.custom_metadata())?.is_some_and(|v| !v.is_empty())
        || checked(schema.features())?.is_some_and(|v| !v.is_empty())
    {
        return Err(err(
            "reopen requires native little-endian IPC without schema extensions",
        ));
    }
    let fields = checked(schema.fields())?.ok_or_else(|| err("missing IPC fields"))?;
    if fields.is_empty() || fields.len() > MAX_COLUMNS {
        return Err(err("IPC must contain 1..64 columns"));
    }
    let mut result = Vec::with_capacity(fields.len());
    let mut names = HashSet::new();
    for field in fields {
        let field = checked(field)?;
        let name = checked(field.name())?.ok_or_else(|| err("missing IPC column name"))?;
        column_name(name)?;
        if !names.insert(name)
            || !checked(field.nullable())?
            || checked(field.dictionary())?.is_some()
            || checked(field.children())?.is_some_and(|v| !v.is_empty())
            || checked(field.custom_metadata())?.is_some_and(|v| !v.is_empty())
        {
            return Err(err(
                "reopen requires native flat IPC fields without dictionaries or extensions",
            ));
        }
        let kind = match checked(field.type_())?.ok_or_else(|| err("missing IPC column dtype"))? {
            ipc::TypeRef::Int(kind)
                if checked(kind.bit_width())? == 64 && checked(kind.is_signed())? =>
            {
                NativeType::Int64
            }
            ipc::TypeRef::FloatingPoint(kind)
                if checked(kind.precision())? == ipc::Precision::Double =>
            {
                NativeType::Float64
            }
            ipc::TypeRef::Bool(_) => NativeType::Boolean,
            ipc::TypeRef::Utf8View(_) => NativeType::StringView,
            ipc::TypeRef::LargeUtf8(_) => NativeType::LargeString,
            _ => {
                return Err(err(
                    "reopen supports only native BioV string/int64/float64/boolean IPC exports",
                ))
            }
        };
        result.push((name.into(), kind));
    }
    Ok(result)
}
fn message(bytes: &[u8], start: usize, metadata_len: usize) -> Result<ipc::MessageRef<'_>> {
    if !(8..=MAX_IPC_METADATA).contains(&metadata_len) {
        return Err(err("IPC message metadata exceeds bounded limit"));
    }
    let data = slice(bytes, start, metadata_len)?;
    if data[..4] != [255; 4] {
        return Err(err("unsupported IPC message framing"));
    }
    let length = nonnegative(i32::from_le_bytes(data[4..8].try_into().unwrap()) as i64)?;
    if length != metadata_len - 8 {
        return Err(err("IPC metadata length prefix does not match block"));
    }
    let message = checked(ipc::MessageRef::read_as_root(&data[8..]))?;
    if checked(message.version())? != ipc::MetadataVersion::V5
        || checked(message.custom_metadata())?.is_some_and(|v| !v.is_empty())
    {
        return Err(err("unsupported IPC message version or metadata"));
    }
    Ok(message)
}
fn preflight_ipc(
    bytes: &[u8],
    declared: &[RecordedColumn],
    rows: usize,
    retained: usize,
) -> Result<()> {
    if bytes.len() < 18 || bytes[..6] != *b"ARROW1" || bytes[bytes.len() - 6..] != *b"ARROW1" {
        return Err(err("invalid native Arrow IPC file signature"));
    }
    let footer_len = nonnegative(i32::from_le_bytes(
        bytes[bytes.len() - 10..bytes.len() - 6].try_into().unwrap(),
    ) as i64)?;
    if footer_len > MAX_IPC_METADATA || footer_len > bytes.len() - 18 {
        return Err(err("IPC footer exceeds bounded metadata limit"));
    }
    let footer_start = bytes.len() - 10 - footer_len;
    let footer = checked(ipc::FooterRef::read_as_root(
        &bytes[footer_start..bytes.len() - 10],
    ))?;
    if checked(footer.version())? != ipc::MetadataVersion::V5
        || checked(footer.dictionaries())?.is_some_and(|v| !v.is_empty())
        || checked(footer.custom_metadata())?.is_some_and(|v| !v.is_empty())
    {
        return Err(err(
            "unsupported IPC footer version, dictionary, or extension",
        ));
    }
    let fields =
        native_schema(checked(footer.schema())?.ok_or_else(|| err("missing IPC schema"))?)?;
    if fields.len() != declared.len()
        || fields.iter().zip(declared).any(|((name, kind), column)| {
            name != &column.name || kind.dtype().to_string() != column.dtype
        })
    {
        return Err(err("IPC schema does not match export record"));
    }
    allocation_budget(rows, fields.len(), retained)?;
    // Also check the initial schema message: it must describe the same supported
    // schema as the footer. Polars uses the footer, other Arrow consumers may not.
    let schema_len =
        nonnegative(i32::from_le_bytes(slice(bytes, 12, 4)?.try_into().unwrap()) as i64)?;
    let schema_meta_len = checked_end(8, schema_len, MAX_IPC_METADATA)?;
    let first = message(bytes, 8, schema_meta_len)?;
    let initial_schema = match checked(first.header())? {
        Some(ipc::MessageHeaderRef::Schema(schema)) => native_schema(schema)?,
        _ => return Err(err("missing initial IPC schema message")),
    };
    if initial_schema != fields || checked(first.body_length())? != 0 {
        return Err(err("IPC initial and footer schemas disagree"));
    }
    let mut previous_end = checked_end(8, schema_meta_len, footer_start)?;
    let batches =
        checked(footer.record_batches())?.ok_or_else(|| err("missing IPC record batches"))?;
    if batches.len() > MAX_BATCHES {
        return Err(err("IPC batch count exceeds bounded limit"));
    }
    let mut total_rows = 0usize;
    let mut total_string_bytes = 0usize;
    let mut decoded_payload = 0usize;
    for block in batches {
        let offset = nonnegative(block.offset())?;
        let metadata_len = nonnegative(i64::from(block.meta_data_length()))?;
        let body_len = nonnegative(block.body_length())?;
        let body_start = checked_end(offset, metadata_len, footer_start)?;
        let end = checked_end(body_start, body_len, footer_start)?;
        if offset < previous_end {
            return Err(err("overlapping or unordered IPC blocks are not supported"));
        }
        previous_end = end;
        let message = message(bytes, offset, metadata_len)?;
        if nonnegative(checked(message.body_length())?)? != body_len {
            return Err(err("IPC body length does not match block"));
        }
        let batch = match checked(message.header())? {
            Some(ipc::MessageHeaderRef::RecordBatch(batch)) => batch,
            _ => return Err(err("expected native IPC record batch")),
        };
        if checked(batch.compression())?.is_some() {
            return Err(err(
                "compressed IPC is not supported for bounded native reopen",
            ));
        }
        let batch_rows = nonnegative(checked(batch.length())?)?;
        total_rows = total_rows
            .checked_add(batch_rows)
            .ok_or_else(|| err("IPC row count overflow"))?;
        if total_rows > rows {
            return Err(err("IPC row count does not match export record"));
        }
        allocation_budget(total_rows, fields.len(), retained)?;
        let nodes = checked(batch.nodes())?.ok_or_else(|| err("missing IPC field nodes"))?;
        if nodes.len() != fields.len() {
            return Err(err("IPC node count does not match schema"));
        }
        let buffers = checked(batch.buffers())?.ok_or_else(|| err("missing IPC buffers"))?;
        if buffers.len() > MAX_BUFFERS {
            return Err(err("IPC buffer count exceeds bounded limit"));
        }
        let counts = checked(batch.variadic_buffer_counts())?;
        let string_fields = fields
            .iter()
            .filter(|(_, kind)| *kind == NativeType::StringView)
            .count();
        if counts.as_ref().map_or(0, |v| v.len()) != string_fields {
            return Err(err("invalid IPC variadic buffer counts"));
        }
        let mut variadic = counts.into_iter().flat_map(|values| values.into_iter());
        let mut buffer_index = 0usize;
        let mut buffer_end = 0usize;
        let mut next_buffer = || -> Result<&[u8]> {
            let buffer = buffers
                .get(buffer_index)
                .ok_or_else(|| err("missing IPC column buffer"))?;
            buffer_index += 1;
            let start = nonnegative(buffer.offset())?;
            let length = nonnegative(buffer.length())?;
            let end = checked_end(start, length, body_len)?;
            if start < buffer_end {
                return Err(err("overlapping IPC buffers are not supported"));
            }
            buffer_end = end;
            slice(bytes, body_start + start, length)
        };
        for ((_, kind), node) in fields.iter().zip(nodes) {
            let length = nonnegative(node.length())?;
            let nulls = nonnegative(node.null_count())?;
            if length != batch_rows || nulls > length {
                return Err(err("invalid IPC field node length or null count"));
            }
            let validity = next_buffer()?;
            let bitmap_bytes = length.div_ceil(8);
            if (validity.len() != bitmap_bytes && !(nulls == 0 && validity.is_empty()))
                || (nulls > 0 && validity.is_empty())
            {
                return Err(err("invalid IPC validity buffer length"));
            }
            if !validity.is_empty() {
                let valid = validity
                    .iter()
                    .enumerate()
                    .map(|(index, byte)| {
                        let used = (length - index * 8).min(8);
                        (byte & ((1u16 << used) - 1) as u8).count_ones() as usize
                    })
                    .sum::<usize>();
                if length - valid != nulls {
                    return Err(err("IPC validity bitmap does not match field null count"));
                }
            }
            // Match memory_charge: StringView estimates count logical string
            // bytes only; its validity is covered by the 17-byte cell charge.
            // Polars drops all-valid bitmaps while reading primitive columns.
            if !kind.is_string() && nulls > 0 {
                decoded_payload = decoded_payload.saturating_add(validity.len());
            }
            let values = next_buffer()?;
            let expected = match kind {
                NativeType::StringView => length * 16,
                NativeType::LargeString => length
                    .checked_add(1)
                    .and_then(|offsets| offsets.checked_mul(8))
                    .ok_or_else(|| err("IPC string offsets allocation overflow"))?,
                NativeType::Int64 | NativeType::Float64 => length * 8,
                NativeType::Boolean => bitmap_bytes,
            };
            if values.len() != expected {
                return Err(err("invalid IPC values buffer length"));
            }
            if !kind.is_string() {
                decoded_payload = decoded_payload.saturating_add(values.len());
            }
            if *kind == NativeType::StringView {
                // Cap logical string bytes as well as physical buffers. A small
                // file can otherwise reference the same huge string many times.
                for view in values.chunks_exact(16) {
                    let length = u32::from_le_bytes(view[..4].try_into().unwrap()) as usize;
                    total_string_bytes = total_string_bytes
                        .checked_add(length)
                        .filter(|v| *v <= MAX_MEMORY_BYTES)
                        .ok_or_else(|| err("IPC string allocation exceeds bounded limit"))?;
                }
                let count = nonnegative(
                    variadic
                        .next()
                        .ok_or_else(|| err("missing IPC variadic count"))?,
                )?;
                if count > MAX_BUFFERS {
                    return Err(err("IPC variadic buffer count exceeds bounded limit"));
                }
                for _ in 0..count {
                    next_buffer()?;
                }
            } else if *kind == NativeType::LargeString {
                // Native compatibility output uses compact i64 offsets and
                // one UTF-8 buffer. Validate before the upstream reader can
                // allocate or convert it to Polars' internal StringView type.
                let data = next_buffer()?;
                let text = checked(std::str::from_utf8(data))?;
                let mut previous = 0usize;
                for (index, offset) in values.chunks_exact(8).enumerate() {
                    let offset = nonnegative(i64::from_le_bytes(offset.try_into().unwrap()))?;
                    if (index == 0 && offset != 0)
                        || offset < previous
                        || offset > data.len()
                        || !text.is_char_boundary(offset)
                    {
                        return Err(err("invalid native IPC string offsets"));
                    }
                    previous = offset;
                }
                if previous != data.len() {
                    return Err(err("IPC string offsets do not cover the values buffer"));
                }
                total_string_bytes = total_string_bytes
                    .checked_add(data.len())
                    .filter(|value| *value <= MAX_MEMORY_BYTES)
                    .ok_or_else(|| err("IPC string allocation exceeds bounded limit"))?;
            }
        }
        if decoded_payload
            .saturating_add(total_string_bytes)
            .saturating_add(rows.saturating_mul(fields.len()).saturating_mul(17))
            .saturating_add(retained)
            > MAX_MEMORY_BYTES
        {
            return Err(err(
                "IPC decoded allocation exceeds retained dataset memory limit",
            ));
        }
        if buffer_index != buffers.len() {
            return Err(err("unexpected extra IPC buffers"));
        }
    }
    if total_rows != rows {
        return Err(err("IPC row count does not match export record"));
    }
    Ok(())
}

#[cfg(test)]
mod adversarial_tests {
    use super::*;
    use std::fs;

    fn fixture() -> (tempfile::TempDir, Value, Vec<u8>) {
        fixture_with_compat(CompatLevel::newest())
    }
    fn fixture_with_compat(compat: CompatLevel) -> (tempfile::TempDir, Value, Vec<u8>) {
        let dir = tempfile::tempdir().unwrap();
        fs::write(dir.path().join("table.csv"), "id,count\n001,1\n002,\n").unwrap();
        let mut store = DatasetStore::new(dir.path(), dir.path()).unwrap();
        let source = store
            .open(
                serde_json::from_value(json!({"path":"table.csv", "schema":{"count":"int64"}}))
                    .unwrap(),
            )
            .unwrap();
        let mut exported = store
            .export(ExportRequest {
                dataset_id: source["dataset_id"].as_str().unwrap().into(),
            })
            .unwrap();
        // Retain explicit old Utf8View fixtures even when production exports
        // use upstream LargeUtf8 compatibility output.
        let mut frame = store
            .get(source["dataset_id"].as_str().unwrap())
            .unwrap()
            .frame
            .clone();
        let mut cursor = Cursor::new(Vec::new());
        IpcWriter::new(&mut cursor)
            .with_compat_level(compat)
            .finish(&mut frame)
            .unwrap();
        let bytes = cursor.into_inner();
        exported["record"]["bytes"] = json!(bytes.len());
        exported["record"]["sha256"] = json!(format!("{:x}", Sha256::digest(&bytes)));
        fs::write(exported["execution_host_path"].as_str().unwrap(), &bytes).unwrap();
        fs::write(
            exported["record_path"].as_str().unwrap(),
            serde_json::to_vec(&exported["record"]).unwrap(),
        )
        .unwrap();
        (dir, exported, bytes)
    }
    fn footer(bytes: &[u8]) -> (usize, ipc::Footer) {
        let len = i32::from_le_bytes(bytes[bytes.len() - 10..bytes.len() - 6].try_into().unwrap())
            as usize;
        let start = bytes.len() - 10 - len;
        (
            start,
            ipc::Footer::try_from(
                ipc::FooterRef::read_as_root(&bytes[start..bytes.len() - 10]).unwrap(),
            )
            .unwrap(),
        )
    }
    fn write_footer(mut prefix: Vec<u8>, footer: &ipc::Footer) -> Vec<u8> {
        let mut builder = ipc::planus::Builder::new();
        let encoded = builder.finish(footer, None);
        prefix.extend_from_slice(encoded);
        prefix.extend_from_slice(&(encoded.len() as i32).to_le_bytes());
        prefix.extend_from_slice(b"ARROW1");
        prefix
    }
    fn mutate_footer(bytes: &[u8], change: impl FnOnce(&mut ipc::Footer)) -> Vec<u8> {
        let (start, mut footer) = footer(bytes);
        change(&mut footer);
        write_footer(bytes[..start].to_vec(), &footer)
    }
    fn mutate_batches(bytes: &[u8], change: impl Fn(&mut ipc::RecordBatch)) -> Vec<u8> {
        let (_, mut footer) = footer(bytes);
        let blocks = footer.record_batches.as_mut().unwrap();
        let mut result = bytes[..blocks[0].offset as usize].to_vec();
        for block in blocks {
            let mut message = ipc::Message::try_from(
                super::message(
                    bytes,
                    block.offset as usize,
                    block.meta_data_length as usize,
                )
                .unwrap(),
            )
            .unwrap();
            let ipc::MessageHeader::RecordBatch(batch) = message.header.as_mut().unwrap() else {
                panic!("expected batch")
            };
            change(batch);
            let mut builder = ipc::planus::Builder::new();
            let encoded = builder.finish(&message, None);
            let length = (encoded.len() + 8).div_ceil(8) * 8;
            let body_start = block.offset as usize + block.meta_data_length as usize;
            block.offset = result.len() as i64;
            block.meta_data_length = length as i32;
            result.extend_from_slice(&[255; 4]);
            result.extend_from_slice(&((length - 8) as i32).to_le_bytes());
            result.extend_from_slice(encoded);
            result.resize(block.offset as usize + length, 0);
            result.extend_from_slice(&bytes[body_start..body_start + block.body_length as usize]);
        }
        result.extend_from_slice(&[255, 255, 255, 255, 0, 0, 0, 0]);
        write_footer(result, &footer)
    }
    fn rewritten_digest_error(dir: &Path, export: &Value, bytes: &[u8]) -> String {
        fs::write(export["execution_host_path"].as_str().unwrap(), bytes).unwrap();
        let mut record = export["record"].clone();
        record["bytes"] = json!(bytes.len());
        record["sha256"] = json!(format!("{:x}", Sha256::digest(bytes)));
        fs::write(
            export["record_path"].as_str().unwrap(),
            serde_json::to_vec(&record).unwrap(),
        )
        .unwrap();
        let filename = Path::new(export["record_path"].as_str().unwrap())
            .file_name()
            .unwrap()
            .to_str()
            .unwrap();
        DatasetStore::new(dir, dir)
            .unwrap()
            .reopen(ReopenRequest {
                record_path: filename.into(),
                preview_rows: 1,
            })
            .unwrap_err()
            .to_string()
    }
    fn first_buffer_range(bytes: &[u8], index: usize) -> std::ops::Range<usize> {
        let (_, footer) = footer(bytes);
        let block = &footer.record_batches.as_ref().unwrap()[0];
        let message = super::message(
            bytes,
            block.offset as usize,
            block.meta_data_length as usize,
        )
        .unwrap();
        let Some(ipc::MessageHeaderRef::RecordBatch(batch)) = message.header().unwrap() else {
            panic!("batch")
        };
        let buffer = batch.buffers().unwrap().unwrap().get(index).unwrap();
        let start =
            block.offset as usize + block.meta_data_length as usize + buffer.offset() as usize;
        start..start + buffer.length() as usize
    }

    #[test]
    fn native_large_utf8_and_prior_string_view_both_reopen() {
        for (compat, is_view) in [
            (CompatLevel::oldest(), false),
            (CompatLevel::newest(), true),
        ] {
            let (dir, export, bytes) = fixture_with_compat(compat);
            let (_, record_footer) = footer(&bytes);
            let physical_type = record_footer.schema.unwrap().fields.unwrap()[0]
                .type_
                .clone()
                .unwrap();
            assert_eq!(matches!(physical_type, ipc::Type::Utf8View(_)), is_view);
            assert_eq!(matches!(physical_type, ipc::Type::LargeUtf8(_)), !is_view);
            let name = Path::new(export["record_path"].as_str().unwrap())
                .file_name()
                .unwrap()
                .to_str()
                .unwrap();
            let mut store = DatasetStore::new(dir.path(), dir.path()).unwrap();
            let opened = store
                .reopen(ReopenRequest {
                    record_path: name.into(),
                    preview_rows: 5,
                })
                .unwrap();
            assert_eq!(
                opened["preview"]["rows"],
                json!([["001", 1], ["002", null]])
            );
            let reexported = store
                .export(ExportRequest {
                    dataset_id: opened["dataset_id"].as_str().unwrap().into(),
                })
                .unwrap();
            let native = fs::read(reexported["execution_host_path"].as_str().unwrap()).unwrap();
            let (_, native_footer) = footer(&native);
            assert!(matches!(
                native_footer.schema.unwrap().fields.unwrap()[0].type_,
                Some(ipc::Type::LargeUtf8(_))
            ));
        }
    }

    #[test]
    fn large_utf8_offsets_and_payload_are_checked_before_decode_even_with_matching_hash() {
        let (dir, export, bytes) = fixture_with_compat(CompatLevel::oldest());
        let offsets = first_buffer_range(&bytes, 1);
        let data = first_buffer_range(&bytes, 2);
        assert_eq!(offsets.len(), 24);
        assert_eq!(data.len(), 6);
        for (offset_values, expected) in [
            ([1, 3, 6], "offsets"), // Native buffers must start at zero.
            ([0, -1, 6], "negative"),
            ([0, 5, 3], "offsets"), // Decreasing offset.
            ([0, 3, 7], "offsets"), // Beyond the physical values buffer.
            ([0, 3, 5], "cover"),   // Unreferenced trailing bytes.
            ([0, i64::MAX, 6], "offsets"),
        ] {
            let mut modified = bytes.clone();
            for (slot, offset) in modified[offsets.clone()]
                .chunks_exact_mut(8)
                .zip(offset_values)
            {
                slot.copy_from_slice(&offset.to_le_bytes());
            }
            let error = rewritten_digest_error(dir.path(), &export, &modified);
            assert!(error.contains(expected), "expected {expected}, got {error}");
        }
        let mut unicode_boundary = bytes.clone();
        unicode_boundary[data.start..data.start + 3].copy_from_slice("基".as_bytes());
        unicode_boundary[offsets.start + 8..offsets.start + 16]
            .copy_from_slice(&1i64.to_le_bytes());
        assert!(rewritten_digest_error(dir.path(), &export, &unicode_boundary).contains("offsets"));
        let mut invalid_utf8 = bytes.clone();
        invalid_utf8[data.start] = 255;
        assert!(rewritten_digest_error(dir.path(), &export, &invalid_utf8).contains("utf-8"));
        for (modified, expected) in [
            (
                mutate_batches(&bytes, |batch| {
                    batch.buffers.as_mut().unwrap()[1].length = 16;
                }),
                "values buffer length",
            ),
            (
                mutate_batches(&bytes, |batch| {
                    batch.buffers.as_mut().unwrap()[2].length = 5;
                }),
                "offsets",
            ),
            (
                mutate_batches(&bytes, |batch| {
                    batch.variadic_buffer_counts = Some(vec![0]);
                }),
                "variadic",
            ),
            (
                mutate_footer(&bytes, |footer| {
                    footer.schema.as_mut().unwrap().fields.as_mut().unwrap()[0].type_ =
                        Some(ipc::Type::Utf8View(Box::default()));
                }),
                "schemas disagree",
            ),
        ] {
            let error = rewritten_digest_error(dir.path(), &export, &modified);
            assert!(error.contains(expected), "expected {expected}, got {error}");
        }
    }
    #[test]
    fn altered_metadata_with_matching_record_hash_is_preflighted_before_decode() {
        let (dir, export, bytes) = fixture();
        let mutations: Vec<(Vec<u8>, &str)> = vec![
            (
                mutate_batches(&bytes, |batch| batch.compression = Some(Box::default())),
                "compressed",
            ),
            (
                mutate_batches(&bytes, |batch| {
                    batch.nodes.as_mut().unwrap()[0].length = i64::MAX
                }),
                "field node",
            ),
            (
                mutate_batches(&bytes, |batch| batch.nodes.as_mut().unwrap()[0].length = -1),
                "negative",
            ),
            (
                mutate_batches(&bytes, |batch| {
                    batch.nodes.as_mut().unwrap()[0].null_count = 99
                }),
                "null count",
            ),
            (
                mutate_batches(&bytes, |batch| {
                    batch.nodes.as_mut().unwrap()[1].null_count = 0
                }),
                "validity bitmap",
            ),
            (
                mutate_batches(&bytes, |batch| batch.length = i64::MAX),
                "row count",
            ),
            (
                mutate_batches(&bytes, |batch| {
                    batch.buffers.as_mut().unwrap()[1].length = i64::MAX
                }),
                "range",
            ),
            (
                mutate_batches(&bytes, |batch| {
                    batch.variadic_buffer_counts.as_mut().unwrap()[0] = i64::MAX
                }),
                "variadic",
            ),
            (
                mutate_footer(&bytes, |footer| {
                    footer.schema.as_mut().unwrap().endianness = ipc::Endianness::Big
                }),
                "little-endian",
            ),
            (
                mutate_footer(&bytes, |footer| {
                    footer.schema.as_mut().unwrap().fields.as_mut().unwrap()[0].nullable = false
                }),
                "native flat",
            ),
            (
                mutate_footer(&bytes, |footer| {
                    let blocks = footer.record_batches.as_mut().unwrap();
                    blocks.push(blocks[0]);
                }),
                "overlapping",
            ),
            (
                mutate_footer(&bytes, |footer| {
                    footer.dictionaries = footer.record_batches.clone()
                }),
                "dictionary",
            ),
        ];
        for (modified, expected) in mutations {
            let error = rewritten_digest_error(dir.path(), &export, &modified);
            assert!(error.contains(expected), "expected {expected}, got {error}");
        }
        let mut bad_footer = bytes.clone();
        let length = bad_footer.len();
        bad_footer[length - 10..length - 6].copy_from_slice(&i32::MAX.to_le_bytes());
        assert!(rewritten_digest_error(dir.path(), &export, &bad_footer).contains("footer"));
        let (_, footer) = footer(&bytes);
        let offset = footer.record_batches.as_ref().unwrap()[0].offset as usize;
        let mut bad_message = bytes.clone();
        bad_message[offset + 4..offset + 8].copy_from_slice(&i32::MAX.to_le_bytes());
        assert!(rewritten_digest_error(dir.path(), &export, &bad_message).contains("length prefix"));
        let block = &footer.record_batches.as_ref().unwrap()[0];
        let message = super::message(&bytes, offset, block.meta_data_length as usize).unwrap();
        let Some(ipc::MessageHeaderRef::RecordBatch(batch)) = message.header().unwrap() else {
            panic!("batch")
        };
        let view_buffer = batch.buffers().unwrap().unwrap().get(1).unwrap();
        let start = offset + block.meta_data_length as usize + view_buffer.offset() as usize;
        let mut inflated = bytes.clone();
        inflated[start..start + 4].copy_from_slice(&((MAX_MEMORY_BYTES + 1) as u32).to_le_bytes());
        assert!(
            rewritten_digest_error(dir.path(), &export, &inflated).contains("string allocation")
        );
    }
    #[test]
    fn preflight_retention_charge_matches_polars_for_nullable_native_columns() {
        let strings = Column::new("id".into(), &[Some("001"), None, Some("")]);
        let unicode = "基因".repeat(64);
        let frames = [
            DataFrame::new(vec![strings.clone()]).unwrap(),
            DataFrame::new(vec![
                Series::new_empty("id".into(), &DataType::String).into()
            ])
            .unwrap(),
            DataFrame::new(vec![Column::new("id".into(), &[None::<&str>, None, None])]).unwrap(),
            DataFrame::new(vec![Column::new(
                "id".into(),
                &[Some(unicode.as_str()), None, Some("🧬")],
            )])
            .unwrap(),
            DataFrame::new(vec![
                strings.clone(),
                Column::new("count".into(), &[Some(1i64), None, Some(-1)]),
                Column::new("ratio".into(), &[Some(1.5f64), None, Some(-1.0)]),
                Column::new("passed".into(), &[Some(true), None, Some(false)]),
            ])
            .unwrap(),
            DataFrame::new(vec![strings.clone()])
                .unwrap()
                .vstack(&DataFrame::new(vec![strings]).unwrap())
                .unwrap(),
        ];
        for frame in frames {
            for compat in [CompatLevel::oldest(), CompatLevel::newest()] {
                let mut frame = frame.clone();
                let mut cursor = Cursor::new(Vec::new());
                IpcWriter::new(&mut cursor)
                    .with_compat_level(compat)
                    .finish(&mut frame)
                    .unwrap();
                let bytes = cursor.into_inner();
                let decoded = IpcReader::new(Cursor::new(bytes.as_slice()))
                    .set_rechunk(false)
                    .finish()
                    .unwrap();
                assert!(frame.equals_missing(&decoded));
                let declared: Vec<RecordedColumn> =
                    serde_json::from_value(schema(&decoded)).unwrap();
                let retained = MAX_MEMORY_BYTES - memory_charge(&decoded);
                preflight_ipc(&bytes, &declared, decoded.height(), retained).unwrap();
                assert!(
                    preflight_ipc(&bytes, &declared, decoded.height(), retained + 1)
                        .unwrap_err()
                        .to_string()
                        .contains("allocation")
                );
            }
        }
    }
    #[test]
    fn native_metadata_reserialization_is_accepted_and_matching_hash_is_not_authenticity() {
        let (dir, export, bytes) = fixture();
        let modified = mutate_batches(&bytes, |_| {});
        fs::write(export["execution_host_path"].as_str().unwrap(), &modified).unwrap();
        let mut record = export["record"].clone();
        record["bytes"] = json!(modified.len());
        record["sha256"] = json!(format!("{:x}", Sha256::digest(&modified)));
        // A coherent record can claim a different historical source. This is
        // explicitly not independently verified merely by matching IPC bytes.
        record["provenance"]["source_sha256"] = json!("a".repeat(64));
        fs::write(
            export["record_path"].as_str().unwrap(),
            serde_json::to_vec(&record).unwrap(),
        )
        .unwrap();
        let filename = Path::new(export["record_path"].as_str().unwrap())
            .file_name()
            .unwrap()
            .to_str()
            .unwrap();
        let opened = DatasetStore::new(dir.path(), dir.path())
            .unwrap()
            .reopen(ReopenRequest {
                record_path: filename.into(),
                preview_rows: 1,
            })
            .unwrap();
        assert_eq!(opened["source_sha256"], "a".repeat(64));
        assert_eq!(
            opened["reopen_verification"]["authenticity"],
            "not_established"
        );
        assert_eq!(opened["reopen_verification"]["original_provenance"], CLAIMS);
    }
}
