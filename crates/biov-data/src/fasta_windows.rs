//! A bounded typed table over one verified native/prepared FASTA sequence.
use super::*;
use biov_prepared::{FastaWindowMetricsRequest, PreparedStore};

const ALGORITHM: &str = "nonoverlapping-fasta-gc-windows";
const REVISION: u32 = 1;
const CANONICAL_POLICY: &str =
    "(G+C)/(A+C+G+T);ascii-case-insensitive;ambiguity-excluded;zero-denominator-null";
const WEIGHTED_POLICY: &str = "biov-core-iupac-dna-gc;all-bases-denominator;N=1/2;B,V=2/3;D,H=1/3";

#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FastaWindowsRequest {
    #[serde(flatten)]
    pub input: FastaWindowMetricsRequest,
    #[serde(default = "default_preview")]
    pub preview_rows: usize,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct SequenceOrigin {
    pub format: String,
    pub reference: String,
    pub snapshot_id: String,
    pub recipe_id: String,
    pub fai_sha256: String,
    pub dictionary_sha256: String,
    pub sequence_id: String,
    pub sequence_length: u64,
    pub window_size: u64,
    pub coordinates: String,
    pub units: String,
    pub canonical_gc_policy: String,
    pub weighted_gc_policy: String,
    pub algorithm: String,
    pub algorithm_revision: u32,
}
impl SequenceOrigin {
    pub(super) fn validate(&self) -> Result<()> {
        let identifier = checked(IdentifierRef::parse(&self.reference))?;
        if identifier.compact_id() != self.reference
            || identifier.namespace() != biov_identifiers::Namespace::RefSeqGcf
            || identifier.version().is_none()
            || self.format != "fasta"
            || self.algorithm != ALGORITHM
            || self.algorithm_revision != REVISION
            || self.coordinates != "0-based-half-open;source-sequence-relative"
            || self.units != "length/counts:bases;gc-fractions:dimensionless"
            || self.canonical_gc_policy != CANONICAL_POLICY
            || self.weighted_gc_policy != WEIGHTED_POLICY
            || self.window_size == 0
            || self.window_size > 1024 * 1024
            || self.sequence_length == 0
            || self.sequence_length > i64::MAX as u64
            || self.sequence_id.is_empty()
            || self.sequence_id.len() > biov_prepared::MAX_ID_BYTES
            || self
                .sequence_id
                .bytes()
                .any(|b| b.is_ascii_whitespace() || !b.is_ascii_graphic())
            || !valid_prefixed_digest(&self.snapshot_id, "sha256-")
            || !valid_prefixed_digest(&self.recipe_id, "sha256-")
            || !valid_digest(&self.fai_sha256)
            || !valid_digest(&self.dictionary_sha256)
        {
            return Err(err("invalid sequence-origin provenance"));
        }
        Ok(())
    }
}
fn valid_digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn valid_prefixed_digest(value: &str, prefix: &str) -> bool {
    value.strip_prefix(prefix).is_some_and(valid_digest)
}

/// Conservative retained row charge before either vectors or Polars allocation.
fn preflight_charge(rows: usize, id_bytes: usize, remaining: usize) -> Result<()> {
    // Nine cells/row charged at the session's 17-byte overhead, plus string
    // payload/offset, five int64s, two float64s, validity and Boolean bytes.
    let row_bytes = 9usize * 17 + id_bytes + 8 + 5 * 8 + 2 * 8 + 9;
    if rows.checked_mul(row_bytes).is_none_or(|n| n > remaining) {
        return Err(err(
            "FASTA window rows/cells/identifiers exceed remaining retained dataset budget",
        ));
    }
    Ok(())
}
impl DatasetStore {
    pub fn fasta_windows(
        &mut self,
        prepared: &PreparedStore,
        request: FastaWindowsRequest,
    ) -> Result<Value> {
        check_rows(request.preview_rows)?;
        if self.datasets.len() >= MAX_DATASETS {
            return Err(err("dataset limit reached; release a dataset first"));
        }
        let selected = checked(prepared.preflight_fasta_window_metrics(request.input))?;
        if selected.length > i64::MAX as u64 {
            return Err(err("sequence length exceeds supported int64 coordinates"));
        }
        let retained = self
            .datasets
            .values()
            .fold(0usize, |n, d| n.saturating_add(memory_charge(&d.frame)));
        let rows = checked(usize::try_from(selected.window_row_count))?;
        preflight_charge(
            rows,
            selected.sequence_id.len(),
            MAX_MEMORY_BYTES.saturating_sub(retained),
        )?;
        let origin = SequenceOrigin {
            format: "fasta".into(),
            reference: selected.reference.clone(),
            snapshot_id: selected.snapshot_id.clone(),
            recipe_id: selected.recipe_id.clone(),
            fai_sha256: selected.fai_sha256.clone(),
            dictionary_sha256: selected.dictionary_sha256.clone(),
            sequence_id: selected.sequence_id.clone(),
            sequence_length: selected.length,
            window_size: selected.window_size,
            coordinates: "0-based-half-open;source-sequence-relative".into(),
            units: "length/counts:bases;gc-fractions:dimensionless".into(),
            canonical_gc_policy: CANONICAL_POLICY.into(),
            weighted_gc_policy: WEIGHTED_POLICY.into(),
            algorithm: ALGORITHM.into(),
            algorithm_revision: REVISION,
        };
        origin.validate()?;
        let source = Path::new(&selected.source_path);
        if selected.source_path.len() > 1024
            || source
                .components()
                .any(|c| !matches!(c, Component::Normal(_)))
        {
            return Err(err(
                "FASTA historical source path exceeds portable table record bounds",
            ));
        }
        let provenance = json!({"source": selected.source_path, "source_bytes": selected.source_bytes, "source_sha256": selected.source_sha256,
            "input_consistency": "verified_native_and_prepared_before_and_after_indexed_read", "operations": [], "declared_schema": {}, "sequence_origin": origin});
        if checked(serde_json::to_vec(&provenance))?.len() > 8 * 1024 {
            return Err(err("provenance exceeds bounded 8 KiB record limit"));
        }
        let metadata = ScientificMetadata {
            identifier: Some(selected.reference.clone()),
            reference: Some(selected.reference.clone()),
            coordinates: Some(origin.coordinates.clone()),
            units: Some(origin.units.clone()),
            species: None,
        };
        validate_metadata(&metadata)?;
        let mut ids = Vec::with_capacity(rows);
        let mut starts = Vec::with_capacity(rows);
        let mut ends = Vec::with_capacity(rows);
        let mut lengths = Vec::with_capacity(rows);
        let mut full = Vec::with_capacity(rows);
        let mut canonical = Vec::with_capacity(rows);
        let mut gc = Vec::with_capacity(rows);
        let mut fractions = Vec::with_capacity(rows);
        let mut weighted = Vec::with_capacity(rows);
        let summary = checked(selected.stream_windows(|row| {
            ids.push(row.sequence_id);
            starts.push(row.start as i64);
            ends.push(row.end as i64);
            lengths.push(row.length as i64);
            full.push(row.is_full_window);
            canonical.push(row.canonical_base_count as i64);
            gc.push(row.gc_base_count as i64);
            fractions.push(row.gc_fraction);
            weighted.push(row.weighted_gc_fraction);
        }))?;
        if ids.len() != rows {
            return Err(err("indexed sequence window count changed"));
        }
        let frame = checked(DataFrame::new(vec![
            Column::new("sequence_id".into(), ids),
            Column::new("start".into(), starts),
            Column::new("end".into(), ends),
            Column::new("length".into(), lengths),
            Column::new("is_full_window".into(), full),
            Column::new("canonical_base_count".into(), canonical),
            Column::new("gc_base_count".into(), gc),
            Column::new("gc_fraction".into(), fractions),
            Column::new("weighted_gc_fraction".into(), weighted),
        ]))?;
        let summary = checked(serde_json::to_value(summary))?;
        let mut result = self.insert(frame, provenance, metadata, None, request.preview_rows)?;
        result["sequence_summary"] = summary;
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn flattened_request_is_strict_and_preserves_opaque_sequence_ids() {
        let request = json!({"reference":"refseq.gcf:GCF_000005845.2", "snapshot_id":format!("sha256-{}", "a".repeat(64)), "source_path":"genome.fna", "recipe_id":format!("sha256-{}", "b".repeat(64)), "sequence_id":"001:chr-a", "window_size":10, "preview_rows":2});
        let parsed: FastaWindowsRequest = serde_json::from_value(request.clone()).unwrap();
        assert_eq!(parsed.input.sequence_id, "001:chr-a");
        assert_eq!(parsed.preview_rows, 2);
        let mut unknown = request;
        unknown["extra"] = json!(true);
        assert!(serde_json::from_value::<FastaWindowsRequest>(unknown).is_err());
    }
    #[test]
    fn window_preflight_includes_identifiers_and_remaining_budget() {
        assert!(preflight_charge(1, 4, 230).is_ok());
        assert!(preflight_charge(1, 4, 229).is_err());
        assert!(preflight_charge(usize::MAX, 4096, MAX_MEMORY_BYTES).is_err());
        assert!(preflight_charge(200_000, 4096, MAX_MEMORY_BYTES).is_err());
    }
}
