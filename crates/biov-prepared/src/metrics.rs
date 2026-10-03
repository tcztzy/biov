//! Bounded indexed metrics over an existing, deeply verified FASTA preparation.
use biov_core::sequence::{nucleotide_gc_counts, Kind, NucleotideGcCounts};
use noodles_core::{Position, Region};
use noodles_fasta::{fai, io::IndexedReader};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{fs, io::BufReader};

use crate::{
    corrupt, fsutil, invalid, io, limit, validate_request, verify_bundle, verify_input_bundle,
    FastaProvenance, Input, PrepareFastaRequest, PreparedError, PreparedStore, IO_BUFFER_BYTES,
    MAX_ID_BYTES, MAX_METADATA_BYTES,
};

/// Maximum bases in any noodles indexed query. No whole sequence is retained.
pub const MAX_WINDOW_BASES: u64 = 1024 * 1024;

/// Analyze one exact sequence from an existing exact prepared FASTA bundle.
/// Window coordinates are zero-based, half-open, with a final short window.
#[derive(Debug, Clone, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FastaWindowMetricsRequest {
    /// Exact canonical versioned native RefSeq reference.
    pub reference: String,
    /// Exact native snapshot digest.
    pub snapshot_id: String,
    /// Source-relative selected genome FASTA file in the native snapshot.
    pub source_path: String,
    /// Exact recipe identity; this operation never creates missing preparations.
    pub recipe_id: String,
    /// Exact opaque FAI name. Case, leading zeroes and punctuation are significant.
    pub sequence_id: String,
    /// Positive window width in bases, at most MAX_WINDOW_BASES.
    pub window_size: u64,
}

impl FastaWindowMetricsRequest {
    pub fn fasta_request(&self) -> PrepareFastaRequest {
        PrepareFastaRequest {
            reference: self.reference.clone(),
            snapshot_id: self.snapshot_id.clone(),
            source_path: self.source_path.clone(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FastaWindowMetric {
    pub sequence_id: String,
    pub start: u64,
    pub end: u64,
    pub length: u64,
    pub is_full_window: bool,
    pub canonical_base_count: u64,
    pub gc_base_count: u64,
    /// Literal G/C divided by canonical A/C/G/T only; null if none exist.
    pub gc_fraction: Option<f64>,
    /// Equal-weight IUPAC GC probability divided by every symbol in the window.
    pub weighted_gc_fraction: f64,
}

/// Whole selected-sequence counts accumulated in the same pass as the windows.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FastaMetricsSummary {
    pub sequence_id: String,
    pub length: u64,
    pub canonical_base_count: u64,
    pub gc_base_count: u64,
    pub gc_fraction: Option<f64>,
    pub weighted_gc_fraction: f64,
}

/// Deeply verified selection with row count available before table allocation.
/// Public fields are descriptive metadata. Streaming uses its private verified
/// request and index record, then revalidates before returning success.
pub struct FastaMetricsSelection {
    pub reference: String,
    pub snapshot_id: String,
    pub recipe_id: String,
    /// Store-relative selected file path, including the exact native snapshot.
    pub source_path: String,
    pub source_sha256: String,
    pub source_bytes: u64,
    pub fai_sha256: String,
    pub dictionary_sha256: String,
    pub sequence_id: String,
    pub length: u64,
    pub window_size: u64,
    pub window_row_count: u64,
    store: PreparedStore,
    request: FastaWindowMetricsRequest,
    input: Input,
    expected: FastaProvenance,
    record: fai::Record,
}

impl PreparedStore {
    /// Verify the native snapshot, independently regenerate and compare the
    /// existing preparation, and select an exact indexed sequence. This reads
    /// only: absent preparations are explicit errors rather than new writes.
    pub fn preflight_fasta_window_metrics(
        &self,
        request: FastaWindowMetricsRequest,
    ) -> Result<FastaMetricsSelection, PreparedError> {
        let fasta_request = request.fasta_request();
        validate_request(&fasta_request)?;
        validate_metrics_request(&request)?;
        fsutil::no_links(&self.root)?;
        let input = self.input(&fasta_request)?;
        let recipe_id = input.recipe.id()?;
        if recipe_id != request.recipe_id {
            return Err(invalid(
                "recipe_id does not match the exact current prepared FASTA recipe",
            ));
        }
        let target = self.root.join("prepared").join(&recipe_id);
        match fs::symlink_metadata(&target) {
            Ok(_) => (),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                return Err(PreparedError::InputUnavailable(
                    "exact prepared FASTA bundle is absent; prepare it explicitly first".into(),
                ));
            }
            Err(error) => return Err(io("inspect prepared metrics target", error)),
        }
        let source = fsutil::open_regular(&self.root.join(&input.relative_path))?;
        let expected = verify_input_bundle(&input, &recipe_id, &target, source)?;
        let fai_identity = expected
            .outputs
            .iter()
            .find(|v| v.path == "sequences.fai")
            .ok_or_else(|| corrupt("verified preparation lacks FAI identity"))?;
        let bytes = fsutil::read_bounded(&target.join("sequences.fai"), MAX_METADATA_BYTES)?;
        if bytes.len() as u64 != fai_identity.bytes
            || format!("{:x}", Sha256::digest(&bytes)) != fai_identity.sha256
        {
            return Err(corrupt("FAI changed after preparation verification"));
        }
        let index = fai::io::Reader::new(bytes.as_slice())
            .read_index()
            .map_err(|_| corrupt("verified FAI cannot be decoded by noodles"))?;
        let record = index
            .as_ref()
            .iter()
            .find(|record| record.name() == request.sequence_id.as_bytes())
            .cloned()
            .ok_or_else(|| invalid("sequence_id is absent from the exact prepared FASTA"))?;
        let length = record.length();
        if length == 0 {
            return Err(corrupt("prepared sequence has zero length"));
        }
        // Explicit region endpoints are usize-backed in noodles. Reject before
        // streaming if the current platform cannot express the selected length.
        usize::try_from(length).map_err(|_| limit("sequence coordinates exceed platform usize"))?;
        let window_row_count =
            length / request.window_size + u64::from(length % request.window_size != 0);
        let fai_sha256 = fai_identity.sha256.clone();
        let dictionary_sha256 = expected
            .outputs
            .iter()
            .find(|v| v.path == "sequences.tsv")
            .ok_or_else(|| corrupt("verified preparation lacks dictionary identity"))?
            .sha256
            .clone();
        self.revalidate(&fasta_request, &input)?;
        Ok(FastaMetricsSelection {
            reference: request.reference.clone(),
            snapshot_id: request.snapshot_id.clone(),
            recipe_id,
            source_path: input.relative_path.clone(),
            source_sha256: input.recipe.source_sha256.clone(),
            source_bytes: input.recipe.source_bytes,
            fai_sha256,
            dictionary_sha256,
            sequence_id: request.sequence_id.clone(),
            length,
            window_size: request.window_size,
            window_row_count,
            store: self.clone(),
            request,
            input,
            expected,
            record,
        })
    }
}

impl FastaMetricsSelection {
    /// Emit bounded windows, never an unbounded whole-contig noodles query. Callers must
    /// discard all callback rows if this returns an error: final native and
    /// prepared verification occurs after emission and before success.
    pub fn stream_windows(
        &self,
        mut emit: impl FnMut(FastaWindowMetric),
    ) -> Result<FastaMetricsSummary, PreparedError> {
        let source = fsutil::open_regular(&self.store.root.join(&self.input.relative_path))?;
        // A single-record index avoids linear lookup through every contig on
        // every window while retaining noodles' own offset and format logic.
        let mut reader = IndexedReader::new(
            BufReader::with_capacity(IO_BUFFER_BYTES, source),
            fai::Index::from(vec![self.record.clone()]),
        );
        let mut total = NucleotideGcCounts::default();
        let mut start = 0u64;
        let length = self.record.length();
        while start < length {
            let window_length = self.request.window_size.min(length - start);
            let end = start + window_length;
            // API coordinates are [start,end); noodles uses [start+1,end].
            // Region::new keeps the identifier opaque rather than parsing ':'
            // or '-' in a textual region expression.
            let first = Position::try_from(
                usize::try_from(start + 1)
                    .map_err(|_| limit("window start exceeds platform coordinates"))?,
            )
            .map_err(|_| corrupt("window start is not a positive noodles position"))?;
            let last = Position::try_from(
                usize::try_from(end)
                    .map_err(|_| limit("window end exceeds platform coordinates"))?,
            )
            .map_err(|_| corrupt("window end is not a positive noodles position"))?;
            let region = Region::new(self.request.sequence_id.as_str(), first..=last);
            let record = reader
                .query(&region)
                .map_err(|e| io("query indexed FASTA window", e))?;
            let sequence: &[u8] = record.sequence().as_ref();
            if sequence.len() as u64 != window_length {
                return Err(corrupt(
                    "indexed FASTA window length differs from verified FAI",
                ));
            }
            let value = std::str::from_utf8(sequence)
                .map_err(|_| corrupt("indexed FASTA window is not ASCII DNA"))?;
            let counts = nucleotide_gc_counts(value, Kind::Dna)
                .map_err(|_| corrupt("indexed FASTA window contains invalid IUPAC DNA"))?;
            total.length += counts.length;
            total.canonical_base_count += counts.canonical_base_count;
            total.gc_base_count += counts.gc_base_count;
            total.weighted_gc_sixths += counts.weighted_gc_sixths;
            emit(FastaWindowMetric {
                sequence_id: self.request.sequence_id.clone(),
                start,
                end,
                length: window_length,
                is_full_window: window_length == self.request.window_size,
                canonical_base_count: counts.canonical_base_count,
                gc_base_count: counts.gc_base_count,
                gc_fraction: counts.gc_fraction(),
                weighted_gc_fraction: counts.weighted_gc_fraction(),
            });
            start = end;
        }
        self.store
            .revalidate(&self.request.fasta_request(), &self.input)?;
        verify_bundle(
            &self
                .store
                .root
                .join("prepared")
                .join(&self.expected.recipe_id),
            &self.expected,
        )?;
        Ok(FastaMetricsSummary {
            sequence_id: self.request.sequence_id.clone(),
            length: total.length,
            canonical_base_count: total.canonical_base_count,
            gc_base_count: total.gc_base_count,
            gc_fraction: total.gc_fraction(),
            weighted_gc_fraction: total.weighted_gc_fraction(),
        })
    }
}

fn validate_metrics_request(request: &FastaWindowMetricsRequest) -> Result<(), PreparedError> {
    let digest = request
        .recipe_id
        .strip_prefix("sha256-")
        .ok_or_else(|| invalid("exact sha256 prepared recipe pin is required"))?;
    if digest.len() != 64
        || !digest
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(invalid(
            "recipe pin must contain 64 lowercase hexadecimal characters",
        ));
    }
    if request.sequence_id.is_empty()
        || request.sequence_id.len() > MAX_ID_BYTES
        || !request
            .sequence_id
            .bytes()
            .all(|b| (b'!'..=b'~').contains(&b))
    {
        return Err(invalid(
            "sequence_id must be a nonempty printable ASCII token of at most 4096 bytes",
        ));
    }
    if request.window_size == 0 {
        return Err(invalid("window_size must be positive"));
    }
    if request.window_size > MAX_WINDOW_BASES {
        return Err(limit("window_size exceeds the 1 MiB indexed-query cap"));
    }
    Ok(())
}
