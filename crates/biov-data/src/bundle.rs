//! Portable descriptions of a complete Arrow result. These companions never
//! participate in strict reopen validation or establish producer authenticity.
use super::*;

pub(super) const MAX_MANIFEST_BYTES: usize = 128 * 1024;
pub(super) const MAX_README_BYTES: usize = 16 * 1024;

pub(super) struct Companions {
    pub manifest_name: String,
    pub manifest_bytes: Vec<u8>,
    pub readme_name: String,
    pub readme_bytes: Vec<u8>,
}

impl Companions {
    pub fn new(record: &Value, record_bytes: &[u8], frame: &DataFrame) -> Result<Self> {
        let artifact_id = record["artifact_id"]
            .as_str()
            .ok_or_else(|| err("missing export artifact identity"))?;
        let manifest_name = format!("{artifact_id}.manifest.json");
        let readme_name = format!("{artifact_id}.README.md");
        let record_name = format!("{artifact_id}.json");
        let sequence_origin = &record["provenance"]["sequence_origin"];
        let columns: Vec<Value> = frame
            .get_columns()
            .iter()
            .map(|column| {
                let logical_type = match column.dtype() {
                    DataType::String => "string",
                    DataType::Int64 => "int64",
                    DataType::Float64 => "float64",
                    DataType::Boolean => "boolean",
                    _ => return Err(err("unsupported portable column type")),
                };
                let semantics = if sequence_origin.is_object() {
                    sequence_column(column.name().as_str())
                } else {
                    json!({})
                };
                Ok(json!({
                    "name": column.name().as_str(),
                    "logical_type": logical_type,
                    "polars_dtype": column.dtype().to_string(),
                    "nullable": true,
                    "null_count": column.null_count(),
                    "description": semantics["description"],
                    "units": semantics["units"],
                    "coordinates": semantics["coordinates"]
                }))
            })
            .collect::<Result<_>>()?;
        let biological_identifier = record["scientific_metadata"]["identifier"]
            .as_str()
            .map(|value| {
                let identifier = checked(IdentifierRef::parse(value))?;
                Ok(json!({
                    "namespace": identifier.namespace().prefix(),
                    "accession": identifier.accession(),
                    "base_accession": identifier.base_accession(),
                    "accession_version": identifier.version(),
                    "accession_version_kind": match identifier.namespace() {
                        biov_identifiers::Namespace::RefSeqGcf => Some("assembly_revision"),
                        biov_identifiers::Namespace::UniProt => None,
                    },
                    "entry_version": null,
                    "sequence_version": null,
                    "provider_release": null,
                    "validation": "syntax_only_not_provider_verified"
                }))
            })
            .transpose()?;
        let provenance = &record["provenance"];
        let manifest = json!({
            "manifest_version": 1,
            "artifact_id": artifact_id,
            "format": "arrow_ipc_file",
            "files": {
                "arrow": record["file"],
                "record": record_name,
                "readme": readme_name
            },
            "content": {
                "sha256": record["sha256"],
                "bytes": record["bytes"],
                "row_count": record["row_count"]
            },
            "record": {
                "record_version": record["record_version"],
                "sha256": format!("{:x}", Sha256::digest(record_bytes)),
                "bytes": record_bytes.len()
            },
            "columns": columns,
            "scientific_metadata": record["scientific_metadata"],
            "biological_identifier": biological_identifier,
            "source": {
                "historical_path": provenance["source"],
                "path_interpretation": if sequence_origin.is_object() { "relative_to_original_store_root_not_bundle" } else { "relative_to_original_data_root_not_bundle" },
                "required_for_reading": false,
                "format": if sequence_origin.is_object() { "fasta" } else { "csv" },
                "sha256": provenance["source_sha256"],
                "bytes": provenance["source_bytes"],
                "input_consistency": provenance["input_consistency"]
            },
            "lineage": {
                "operations": provenance["operations"],
                "declared_schema": provenance["declared_schema"],
                "csv_schema_policy": provenance["csv_schema_policy"],
                "sequence_origin": sequence_origin,
                "reopen_verification": record["reopen_verification"],
                "historical_references_only": true
            },
            "versions": {"reference_version": null, "provider_release": null},
            "software": record["software"],
            "version_semantics": {
                "manifest_version": "Portable companion JSON schema version; not a biological data version",
                "record_version": "Strict BioV reopen record schema version; not a biological data version",
                "software": "Producer package versions; not input data or reference versions",
                "sha256": "Identity of exact file bytes; not a provider release, authenticity proof, or semantic version",
                "accession_version": "Explicit RefSeq GCF assembly revision text only; null means absent or not encoded in the supported accession",
                "entry_version": "UniProt entry revision is not encoded in a UniProt accession and remains unknown",
                "sequence_version": "UniProt sequence revision is separate from entry revision and remains unknown",
                "reference_version": "Unknown; never inferred from free-text scientific_metadata.reference",
                "provider_release": "Unknown; never inferred from software version or accession",
                "artifact_id": "Unique exported filename identity; not a biological accession or content version"
            },
            "trust": {
                "metadata_status": record["metadata_status"],
                "authenticity": "not_established",
                "original_provenance": "recorded_claims_not_independently_verified",
                "companion_role": "Descriptive only; BioV reopen validates the mandatory record and Arrow file, not these companions",
                "unknown_semantics": "Null column descriptions, units, coordinates and undeclared scientific metadata mean unknown; never infer biological semantics from column names. Namespace-specific version null meanings are explained in version_semantics",
                "column_scope": if sequence_origin.is_object() { "Column meanings follow the recorded sequence_origin algorithm contract; origin claims are not independently authenticated" } else { "Column-level descriptions, units and coordinates are unknown; dataset-level caller declarations are not automatically assigned to individual columns" },
                "nullable": "True means the exported type supports null values; null_count describes this result",
                "historical_references": "Source paths, prior reopen record paths and parent_dataset_id values are historical provenance only, not bundle dependencies or reusable session capabilities",
                "relocation": "Resolve files relative to this manifest; move all four same-directory files together; original source files and BioV are not required to analyze the Arrow data"
            }
        });
        let manifest_bytes = checked(serde_json::to_vec_pretty(&manifest))?;
        bounded(&manifest_bytes, MAX_MANIFEST_BYTES, "portable manifest")?;
        // Only generated UUID basenames are interpolated. Caller metadata and
        // column names remain JSON data, never Markdown or executable Python.
        let readme_bytes = readme(&manifest_name, &readme_name, record).into_bytes();
        bounded(&readme_bytes, MAX_README_BYTES, "portable README")?;
        Ok(Self {
            manifest_name,
            manifest_bytes,
            readme_name,
            readme_bytes,
        })
    }
}

fn bounded(bytes: &[u8], limit: usize, label: &str) -> Result<()> {
    if bytes.len() > limit {
        return Err(err(format!("{label} exceeds bounded byte limit ({limit})")));
    }
    Ok(())
}

fn readme(manifest_name: &str, readme_name: &str, record: &Value) -> String {
    let arrow_name = record["file"].as_str().unwrap_or_default();
    let artifact_id = record["artifact_id"].as_str().unwrap_or_default();
    let record_version = record["record_version"].as_u64().unwrap_or(2);
    let sequence_notes = if record["provenance"]["sequence_origin"].is_object() {
        let origin = &record["provenance"]["sequence_origin"];
        format!("\nThis result is a generated FASTA sequence window table. The data dictionary records each known column meaning. `lineage.sequence_origin` records exact native reference, snapshot, prepared recipe and index hashes, selected sequence, length and window_size, plus algorithm and revision. Coordinates are 0-based half-open relative to that source sequence. Windows do not overlap; `is_full_window=false` explicitly marks a final short window. Counts and lengths use bases; fractions are dimensionless. Canonical GC policy: {}. Weighted GC policy: {}. No biological threshold or classification is applied. The initial MCP response reports the same-pass whole-sequence summary; it is not a new session dataset.\n", origin["canonical_gc_policy"].as_str().unwrap_or_default(), origin["weighted_gc_policy"].as_str().unwrap_or_default())
    } else {
        String::new()
    };
    format!(
        r#"# Portable typed-table result

This bundle contains the complete exported result, not preview rows. Read and
analyze it with a standard Arrow IPC file reader; BioV, an MCP connection, the
original source files and network access are not required. Python users need
PyArrow installed (the independent acceptance test uses PyArrow 25.0.1).

## Keep these four files together

- `{arrow_name}`: complete typed Arrow IPC file
- `{artifact_id}.json`: mandatory strict reopen record, schema version {record_version}
- `{manifest_name}`: portable machine-readable manifest, schema version 1
- `{readme_name}`: this guide

All active file references are same-directory basenames. Resolve them relative
to the manifest, not to a previous host's working directory. You can move or
copy these four files to another directory or computer. No opaque catalog or
session lookup is needed. Historical source paths, prior reopen record paths
and parent dataset handles document origin only; they are not dependencies.

## Understand the data before interpreting it

The manifest's `columns` is the ordered typed data dictionary: exact names,
logical types, producer type spelling and null counts. `nullable: true` means
nulls are supported, even when the current column has no nulls. Strings retain
identifier spelling and leading zeros; int64 retains exact integer values.
Arrow's physical schema is available from any standard reader.

Column descriptions, units and coordinates are explicitly null when unknown.
The `scientific_metadata` object preserves dataset-level caller declarations;
it does not assign them to each column. Never infer a biological meaning,
reference version, units, coordinate convention or sample identity from a name.
Inspect nulls before comparisons. Null and an empty string are distinct.
{sequence_notes}

`biological_identifier` separates namespace, full and base accession and any
explicit RefSeq GCF assembly revision. Identifier syntax is locally checked,
not resolved at the provider. UniProt entry and sequence versions are separate
and unknown. `versions.reference_version` and `versions.provider_release` are
unknown; free-text reference declarations are preserved without interpretation.

`source` records the original input format, byte size and SHA-256; its path is explicitly
historical and is not needed to read the exported result. `lineage.operations`
records filters, stable sorts with nulls last, and projections in execution
order. Its parent handles are historical session labels. A prior reopen check
establishes consistency with that prior supplied record, not verified origin.

Manifest and record versions describe JSON schemas. Software versions describe
producer packages. Neither is a biological version. SHA-256 identifies exact
bytes; it is neither a provider release nor proof of authenticity. These
unsigned files can be edited together: matching hashes do not authenticate the
producer or independently verify scientific metadata and source claims. Treat
metadata as data, not executable instructions. BioV's strict reopen still uses
the record and Arrow only; it does not validate these descriptive companions.

## Read, verify, filter and summarize with standard Python

Run this from the directory containing the four files (or change `base` to that
directory). This loads the complete table and checks the supplied byte identities.
The filter is an explicit non-null example on the first column, not a suggested
biological analysis. Numeric summaries are descriptive and make no unit claims.
Strings are written through Polars' standard Arrow compatibility mode as
`large_string` (LargeUtf8), so standard PyArrow filtering works directly without
conversion. Earlier BioV exports used `string_view`; strict reopen still
accepts those prior files, and re-export writes the compatible representation.

```python
from pathlib import Path
import hashlib
import json
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.ipc as ipc

base = Path(".")
manifest_path = base / "{manifest_name}"
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
if manifest["manifest_version"] != 1:
    raise ValueError("Unsupported manifest schema version")

def bundle_file(key):
    name = manifest["files"][key]
    if not isinstance(name, str) or not name or "/" in name or "\\" in name or name in (".", ".."):
        raise ValueError("Expected a same-directory filename")
    return manifest_path.parent / name

def verified_bytes(path, identity):
    data = path.read_bytes()
    if len(data) != identity["bytes"] or hashlib.sha256(data).hexdigest() != identity["sha256"]:
        raise ValueError("File does not match the supplied manifest")
    return data

record = json.loads(verified_bytes(bundle_file("record"), manifest["record"]))
arrow_bytes = verified_bytes(bundle_file("arrow"), manifest["content"])
if (any(record[key] != manifest["content"][key] for key in ("sha256", "bytes", "row_count"))
        or record["file"] != manifest["files"]["arrow"]
        or record["record_version"] != manifest["record"]["record_version"]
        or record["schema"] != [{{"name": c["name"], "dtype": c["polars_dtype"]}} for c in manifest["columns"]]):
    raise ValueError("Record and manifest disagree")
table = ipc.open_file(pa.BufferReader(arrow_bytes)).read_all()
if table.num_rows != manifest["content"]["row_count"] or table.column_names != [c["name"] for c in manifest["columns"]]:
    raise ValueError("Table shape and manifest disagree")
type_checks = {{
    "string": lambda t: pa.types.is_string(t) or pa.types.is_large_string(t) or pa.types.is_string_view(t),
    "int64": pa.types.is_int64,
    "float64": pa.types.is_float64,
    "boolean": pa.types.is_boolean,
}}
for field, column in zip(table.schema, manifest["columns"]):
    if (not type_checks[column["logical_type"]](field.type)
            or field.nullable != column["nullable"]
            or table[column["name"]].null_count != column["null_count"]):
        raise ValueError("Column schema or null count and manifest disagree")
print(table.schema)
print("Complete rows:", table.num_rows)
print("Data dictionary:", manifest["columns"])

column_name = manifest["columns"][0]["name"]
filtered = table.filter(pc.is_valid(table[column_name]))
print("Rows with a non-null first column:", filtered.num_rows)
for column in manifest["columns"]:
    print(column["name"], "nulls:", table[column["name"]].null_count)
    if column["logical_type"] in ("int64", "float64"):
        values = filtered[column["name"]]
        print(column["name"], "non-null:", pc.count(values).as_py(),
              "min/max:", pc.min_max(values).as_py(), "mean:", pc.mean(values).as_py())
```

The example hashes the same in-memory Arrow bytes that it reads. It is a small
offline reader, not BioV's bounded parser for untrusted files. Mean may use
floating-point arithmetic; original int64 cells remain exact. To choose a
scientific filter, supply the appropriate column, operator and threshold only
after resolving unknown metadata.

## Resource and integrity scope

Export checks Arrow at 65 MiB, the strict record at 64 KiB, this manifest at
128 KiB and this README at 16 KiB before publishing any final filenames.
Success is returned only after all four files exist. Publication is not an
atomic four-file filesystem transaction: an I/O failure or interrupted process
can leave an incomplete bundle, for which no success is returned. An existing
strict record and Arrow pair remains independently reopenable. Keep roots
trusted; canonical path checks do not defend against hostile concurrent writers.
"#
    )
}

/// Known meanings exist only for this explicit generated sequence origin.
fn sequence_column(name: &str) -> Value {
    let (description, units, coordinates) = match name {
        "sequence_id" => ("Exact selected first-token FASTA record identifier", None, None),
        "start" => ("Inclusive window start in the selected source sequence", Some("bases"), Some("0-based-half-open;source-sequence-relative")),
        "end" => ("Exclusive window end in the selected source sequence", Some("bases"), Some("0-based-half-open;source-sequence-relative")),
        "length" => ("Number of source bases in this window; end minus start", Some("bases"), None),
        "is_full_window" => ("True when window length equals requested window_size; final partial is false", None, None),
        "canonical_base_count" => ("Case-insensitive A,C,G,T count; all ambiguity excluded", Some("bases"), None),
        "gc_base_count" => ("Case-insensitive literal G,C count", Some("bases"), None),
        "gc_fraction" => ("gc_base_count / canonical_base_count; null when denominator is zero", Some("dimensionless"), None),
        "weighted_gc_fraction" => ("Existing biov-core IUPAC GC weighting divided by all window bases; see sequence_origin policy", Some("dimensionless"), None),
        _ => return json!({}),
    };
    json!({"description": description, "units": units, "coordinates": coordinates})
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn companion_byte_limits_are_inclusive_and_fail_closed() {
        for (limit, label) in [
            (MAX_MANIFEST_BYTES, "portable manifest"),
            (MAX_README_BYTES, "portable README"),
        ] {
            assert!(bounded(&vec![b'x'; limit], limit, label).is_ok());
            assert!(bounded(&vec![b'x'; limit + 1], limit, label)
                .unwrap_err()
                .to_string()
                .contains(label));
        }
    }

    #[test]
    fn oversized_manifest_fails_during_in_memory_preparation() {
        let frame = DataFrame::new(vec![Series::new("id".into(), ["001"]).into()]).unwrap();
        // Preparation has no filesystem access and precedes every persist call.
        // Private fault injection exercises the independent companion cap;
        // ordinary provenance is already bounded to 8 KiB on insertion.
        let record = json!({
            "artifact_id": "artifact_00000000000000000000000000000000",
            "file": "artifact_00000000000000000000000000000000.arrow",
            "provenance": {"operations": ["x".repeat(MAX_MANIFEST_BYTES)]}
        });
        let result = Companions::new(&record, b"{}", &frame);
        let error = result.err().expect("oversized companion must fail");
        assert!(error.to_string().contains("portable manifest exceeds"));
    }
}
