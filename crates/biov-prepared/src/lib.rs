//! Immutable, portable FASTA indices over exact verified native snapshots.
//!
//! The runtime-independent library keeps native source bytes untouched and uses
//! pinned noodles-fasta for format decoding and conventional FAI serialization.
//! See the crate README for limits, recipe identity and portable closure.
mod fasta;
mod fsutil;
mod types;
pub use types::*;

use biov_identifiers::{IdentifierRef, Namespace};
use biov_storage::{
    EntryKind, NativeStore, Receipt, Resolution, ResolveRequest, ResolvedFile, SnapshotSummary,
};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::BufWriter,
    path::{Path, PathBuf},
};

const RECIPE_PREFIX: &[u8] = b"biov-prepared-fasta-recipe-v1\0";

/// A handle to an existing trusted native store. Construction is read-only.
/// Path checks do not sandbox hostile concurrent local filesystem writers.
#[derive(Debug, Clone)]
pub struct PreparedStore {
    root: PathBuf,
    native: NativeStore,
}

struct Input {
    recipe: FastaRecipe,
    snapshot: SnapshotSummary,
    relative_path: String,
}

impl PreparedStore {
    pub fn new(store_root: impl AsRef<Path>) -> Result<Self, PreparedError> {
        let root = fsutil::root(store_root.as_ref())?;
        let native = NativeStore::new(&root)?;
        Ok(Self { root, native })
    }

    /// Prepare or deeply verify reuse of one already registered plain genomic
    /// FASTA. An exact biological version, snapshot pin and selected file are
    /// mandatory. This performs full native verification and sequential indexing,
    /// including on reuse; it does not promise cheap metadata-only lookup.
    pub fn prepare_fasta(
        &self,
        request: PrepareFastaRequest,
    ) -> Result<PreparedFasta, PreparedError> {
        validate_request(&request)?;
        fsutil::no_links(&self.root)?;
        let input = self.input(&request)?;
        let recipe_id = input.recipe.id()?;
        let _lock = fsutil::lock(&self.root, &recipe_id)?;
        let parent = fsutil::ensure_dir(&self.root, "prepared")?;
        let target = parent.join(&recipe_id);
        let exists = match fs::symlink_metadata(&target) {
            Ok(_) => true,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => false,
            Err(error) => return Err(io("inspect prepared target", error)),
        };
        let source = fsutil::open_regular(&self.root.join(&input.relative_path))?;
        if exists {
            // Regenerating into hash-only sinks verifies lineage as well as saved
            // hashes: editing both an index and its provenance cannot invent a hit.
            let indexed = fasta::index(source, std::io::sink(), std::io::sink())?;
            let expected = provenance(&input, &recipe_id, indexed)?;
            verify_bundle(&target, &expected)?;
            self.revalidate(&request, &input)?;
            return self.result(&input, &expected, true);
        }
        let stages = fsutil::ensure_dir(&self.root, ".staging")?;
        let stage = tempfile::Builder::new()
            .prefix("prepare-fasta-")
            .tempdir_in(&stages)
            .map_err(|e| io("create prepared staging directory", e))?;
        let mut fai = BufWriter::with_capacity(
            IO_BUFFER_BYTES,
            fsutil::create_file(&stage.path().join("sequences.fai"))?,
        );
        let mut dictionary = BufWriter::with_capacity(
            IO_BUFFER_BYTES,
            fsutil::create_file(&stage.path().join("sequences.tsv"))?,
        );
        let indexed = fasta::index(source, &mut fai, &mut dictionary)?;
        fai.get_ref()
            .sync_all()
            .map_err(|e| io("sync FASTA index", e))?;
        dictionary
            .get_ref()
            .sync_all()
            .map_err(|e| io("sync sequence dictionary", e))?;
        drop((fai, dictionary));
        let expected = provenance(&input, &recipe_id, indexed)?;
        let mut result = self.result(&input, &expected, false)?;
        let readme = prepared_readme(&expected);
        fsutil::write_new(&stage.path().join("README.md"), readme.as_bytes())?;
        let json = serde_json::to_vec_pretty(&expected)
            .map_err(|_| limit("cannot serialize provenance"))?;
        if json.len() > MAX_PROVENANCE_BYTES {
            return Err(limit("provenance exceeds 64 KiB"));
        }
        fsutil::write_new(&stage.path().join("provenance.json"), &json)?;
        verify_bundle(stage.path(), &expected)?;
        // NativeStore rehashes the complete pinned tree. Combined with the hash
        // of the very stream indexed above, this detects ordinary in-flight edits.
        self.revalidate(&request, &input)?;
        fsutil::sync_dir(stage.path())?;
        fsutil::no_links(&parent)?;
        let reused = match fsutil::publish(stage.path(), &target) {
            Ok(()) => false,
            Err(PreparedError::Io {
                kind: std::io::ErrorKind::AlreadyExists,
                ..
            }) => {
                verify_bundle(&target, &expected)?;
                true
            }
            Err(error) => return Err(error),
        };
        fsutil::sync_dir(&parent)?;
        fsutil::sync_dir(&stages)?;
        result.reused = reused;
        Ok(result)
    }

    fn input(&self, request: &PrepareFastaRequest) -> Result<Input, PreparedError> {
        let resolution = self.native.resolve(ResolveRequest {
            reference: request.reference.clone(),
            representation: "genome_fasta".into(),
            snapshot_id: Some(request.snapshot_id.clone()),
            scope: Some("assembly".into()),
        })?;
        let (snapshot, paths) = match resolution {
            Resolution::Ready { snapshot, paths } => (snapshot, paths),
            Resolution::Corrupt { detail, .. } => return Err(corrupt(&detail)),
            Resolution::Miss => {
                return Err(PreparedError::InputUnavailable(
                    "exact native snapshot is absent".into(),
                ))
            }
            Resolution::Unavailable { .. } => {
                return Err(PreparedError::InputUnavailable(
                    "genome_fasta is absent from the native package".into(),
                ))
            }
            Resolution::Ambiguous { .. } => {
                return Err(PreparedError::InputUnavailable(
                    "native selection is ambiguous".into(),
                ))
            }
        };
        if snapshot.canonical_ref != request.reference
            || snapshot.snapshot_id != request.snapshot_id
            || snapshot.scope != "assembly"
        {
            return Err(corrupt(
                "native snapshot resolution disagrees with exact request",
            ));
        }
        let relative_path = format!("{}/source/{}", snapshot.snapshot_path, request.source_path);
        if !paths.iter().any(|p| p.relative_path == relative_path) {
            return Err(invalid(
                "source_path is not a catalog-selected genome_fasta file in the exact snapshot",
            ));
        }
        let receipt: Receipt = serde_json::from_slice(&fsutil::read_bounded(
            &self.root.join(&snapshot.receipt_path),
            biov_storage::MAX_METADATA_BYTES,
        )?)
        .map_err(|_| corrupt("malformed native receipt"))?;
        if receipt.snapshot_id != snapshot.snapshot_id
            || receipt.source_content_sha256 != snapshot.source_content_sha256
        {
            return Err(corrupt("native receipt changed after resolution"));
        }
        let entry = receipt
            .inventory
            .iter()
            .find(|entry| entry.path == request.source_path && entry.kind == EntryKind::File)
            .ok_or_else(|| corrupt("selected native file is absent from verified inventory"))?;
        let source_sha256 = entry
            .sha256
            .clone()
            .ok_or_else(|| corrupt("selected native file lacks a hash"))?;
        let recipe = FastaRecipe {
            schema_version: 1, algorithm: "plain-genomic-fasta-fai".into(), algorithm_revision: 1,
            implementation: "biov-prepared".into(), implementation_version: env!("CARGO_PKG_VERSION").into(),
            library: "noodles-fasta".into(), library_version: "0.66.0".into(),
            reference: request.reference.clone(), snapshot_id: request.snapshot_id.clone(),
            source_content_sha256: snapshot.source_content_sha256.clone(), source_path: request.source_path.clone(),
            source_sha256, source_bytes: entry.bytes,
            parameters: FastaParameters {
                input_format: "uncompressed-fasta-lf-or-crlf".into(),
                identifier_convention: "first-ascii-whitespace-delimited-token-after->;exact;unique".into(),
                sequence_byte_policy: "iupac-dna-ascii-case-preserved;no-blank-lines;printable-ascii-or-tab-definitions".into(),
                dictionary_schema: "tsv-v1:sequence_id(string),length(u64-bases);source-order".into(),
                max_physical_line_bytes: MAX_LINE_BYTES, max_identifier_bytes: MAX_ID_BYTES,
                max_records: MAX_RECORDS, max_metadata_bytes: MAX_METADATA_BYTES,
            },
        };
        Ok(Input {
            recipe,
            snapshot,
            relative_path,
        })
    }

    fn revalidate(
        &self,
        request: &PrepareFastaRequest,
        expected: &Input,
    ) -> Result<(), PreparedError> {
        let actual = self.input(request)?;
        if actual.recipe != expected.recipe || actual.relative_path != expected.relative_path {
            return Err(corrupt("native source changed during preparation"));
        }
        Ok(())
    }

    fn result(
        &self,
        input: &Input,
        provenance: &FastaProvenance,
        reused: bool,
    ) -> Result<PreparedFasta, PreparedError> {
        let base = format!("prepared/{}", provenance.recipe_id);
        let result = PreparedFasta {
            recipe_id: provenance.recipe_id.clone(),
            reused,
            reference: input.recipe.reference.clone(),
            snapshot_id: input.recipe.snapshot_id.clone(),
            sequence_count: provenance.sequence_count,
            total_bases: provenance.total_bases,
            fasta: self.path(&input.relative_path)?,
            fai: self.path(&format!("{base}/sequences.fai"))?,
            dictionary: self.path(&format!("{base}/sequences.tsv"))?,
            provenance: self.path(&format!("{base}/provenance.json"))?,
            readme: self.path(&format!("{base}/README.md"))?,
        };
        bounded_response(result)
    }
    fn path(&self, relative: &str) -> Result<ResolvedFile, PreparedError> {
        Ok(ResolvedFile {
            relative_path: relative.into(),
            execution_host_path: self
                .root
                .join(relative)
                .to_str()
                .ok_or_else(|| invalid("execution host path is not UTF-8"))?
                .into(),
        })
    }
}

impl FastaRecipe {
    /// `sha256-` of the documented prefix followed by compact UTF-8 JSON in this
    /// struct's declared field order. No output hashes or host paths are included.
    pub fn id(&self) -> Result<String, PreparedError> {
        let bytes =
            serde_json::to_vec(self).map_err(|_| invalid("cannot serialize FASTA recipe"))?;
        if bytes.len() > MAX_PROVENANCE_BYTES {
            return Err(limit("recipe exceeds 64 KiB"));
        }
        let mut hash = Sha256::new();
        hash.update(RECIPE_PREFIX);
        hash.update(bytes);
        Ok(format!("sha256-{:x}", hash.finalize()))
    }
}

fn validate_request(request: &PrepareFastaRequest) -> Result<(), PreparedError> {
    if request.reference.len() > 256 {
        return Err(invalid("canonical reference exceeds 256 bytes"));
    }
    let reference = IdentifierRef::parse(&request.reference)
        .map_err(|_| invalid("expected canonical versioned refseq.gcf reference"))?;
    if reference.namespace() != Namespace::RefSeqGcf
        || reference.version().is_none()
        || reference.compact_id() != request.reference
    {
        return Err(invalid(
            "expected exact canonical versioned refseq.gcf:GCF_ reference",
        ));
    }
    let digest = request
        .snapshot_id
        .strip_prefix("sha256-")
        .ok_or_else(|| invalid("exact sha256 snapshot pin is required"))?;
    if digest.len() != 64
        || !digest
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(invalid(
            "snapshot pin must contain 64 lowercase hexadecimal characters",
        ));
    }
    fsutil::relative(&request.source_path)?;
    let lower = request.source_path.to_ascii_lowercase();
    if [".gz", ".bgz", ".bgzf", ".bz2", ".xz", ".zip", ".zst"]
        .iter()
        .any(|extension| lower.ends_with(extension))
    {
        return Err(PreparedError::InvalidFasta(
            "compressed FASTA is not supported by this recipe".into(),
        ));
    }
    Ok(())
}

fn provenance(
    input: &Input,
    recipe_id: &str,
    indexed: fasta::Indexed,
) -> Result<FastaProvenance, PreparedError> {
    if indexed.source_sha256 != input.recipe.source_sha256
        || indexed.source_bytes != input.recipe.source_bytes
    {
        return Err(corrupt(
            "bytes consumed by the indexer differ from the verified source inventory",
        ));
    }
    let mut result = FastaProvenance {
        schema_version: 1,
        recipe_id: recipe_id.into(),
        recipe: input.recipe.clone(),
        source_relative_path: format!("../../{}", input.relative_path),
        source_snapshot_relative_path: format!("../../{}", input.snapshot.snapshot_path),
        sequence_count: indexed.sequence_count,
        total_bases: indexed.total_bases,
        outputs: vec![indexed.fai, indexed.dictionary],
    };
    let readme = prepared_readme(&result);
    if readme.len() > MAX_PROVENANCE_BYTES {
        return Err(limit("README exceeds 64 KiB"));
    }
    result.outputs.push(OutputIdentity {
        path: "README.md".into(),
        bytes: readme.len() as u64,
        sha256: format!("{:x}", Sha256::digest(readme.as_bytes())),
    });
    Ok(result)
}

fn verify_bundle(path: &Path, expected: &FastaProvenance) -> Result<(), PreparedError> {
    fsutil::verify_members(path)?;
    let actual: FastaProvenance = serde_json::from_slice(&fsutil::read_bounded(
        &path.join("provenance.json"),
        MAX_PROVENANCE_BYTES,
    )?)
    .map_err(|_| corrupt("malformed or unsupported prepared provenance"))?;
    if &actual != expected {
        return Err(corrupt(
            "prepared provenance disagrees with exact recipe or independently regenerated outputs",
        ));
    }
    for output in &expected.outputs {
        if fsutil::hash_file(&path.join(&output.path), &output.path, MAX_METADATA_BYTES)? != *output
        {
            return Err(corrupt("prepared output checksum or length mismatch"));
        }
    }
    Ok(())
}

fn bounded_response<T: Serialize>(result: T) -> Result<T, PreparedError> {
    let bytes =
        serde_json::to_vec(&result).map_err(|_| limit("cannot serialize prepared response"))?;
    if bytes.len() > MAX_PROVENANCE_BYTES {
        return Err(limit("prepared response exceeds 64 KiB"));
    }
    Ok(result)
}

fn prepared_readme(p: &FastaProvenance) -> String {
    format!(
        r#"# Prepared genomic FASTA index

Reference: {reference}
Native snapshot: {snapshot}
Recipe: {recipe}

This directory contains standard `sequences.fai` and `sequences.tsv`, produced by
biov-prepared {version}, algorithm plain-genomic-fasta-fai revision 1, using
noodles-fasta 0.66.0. Native FASTA bytes are unchanged. See provenance.json for
exact relative paths, byte sizes, SHA-256 hashes, recipe inputs and fixed limits.

## Portable closure and known meaning

The source FASTA is at the `source_relative_path` in provenance.json. It belongs
to the complete native snapshot at `source_snapshot_relative_path`. Copy that
ENTIRE snapshot (native source/, acquisition.json, checksums.sha256 and README.md)
and this ENTIRE prepared directory, preserving their artifacts/... and
prepared/... layout beneath any new root. Copying this prepared directory alone
is incomplete. No original import path, catalog database, BioV install or network
is needed. This is a local recipe, not a general export/materialization command.
The store's native registration README describes its native component; this
separate prepared/ layer may coexist with it and never rewrites native wrappers.

`sequences.fai` is the conventional five-column FASTA index: sequence name,
sequence length in bases, zero-based byte offset of its first base, bases per
full sequence line, and bytes per full sequence line including the terminator.
Identifiers are the exact first whitespace-delimited token immediately after >;
descriptions remain only in the source FASTA. Order is source record order.
`sequences.tsv` is UTF-8 tab-separated text with a header and two columns:
sequence_id (opaque case-preserved string) and length (unsigned integer bases).
Leading zeroes and punctuation in identifiers are significant. Lengths are
literal base counts; this table contains no genomic feature coordinates.
Assembly identity comes from the native RefSeq catalog validation. Chromosome
roles, annotation version, taxonomy, sample, units beyond base counts, acquisition
time and scientific QC are not invented here; consult native metadata for facts.
Hashes establish consistency, not producer authenticity or scientific validity.

## Supported input and operational limits

Plain uncompressed FASTA only, with LF or CRLF physical terminators (final LF
optional), nonempty records and no blank lines. Sequence bytes must be IUPAC DNA
letters A C G T R Y S W K M B D H V N in either case, preserved without normalization.
Definition lines are printable ASCII, with horizontal tabs permitted. Bare CR,
non-ASCII data, gaps, sequence whitespace and duplicate IDs are rejected. Upstream
noodles enforces FAI-compatible wrapping: all nonfinal sequence lines within a
record have equal byte width/base count; the final line may be shorter.
Maximum physical line: 1 MiB including newline; maximum identifier: 4096 bytes;
maximum records: 100000. Each FAI/TSV and the charged retained identifier set has
a 16 MiB limit. Provenance and this README have 64 KiB limits. These independent
bounds are not an aggregate peak-memory guarantee or a total sequence-size cap.
Indexing uses a 64 KiB input buffer, one capped physical line and bounded metadata.
Reuse rehashes the native snapshot, regenerates index/dictionary hashes from the
same indexed input stream, verifies every output and revalidates the source.
Trusted local roots are not a sandbox against hostile concurrent filesystem edits.
Only tested Linux behavior is validated; no distributed filesystem claim is made.

## Ordinary reader example

Install pysam separately in a clean reader environment. From this directory, this
reads every complete record, filters for length >= 4 and summarizes literal G/C.
The explicit external index path is essential: no index is written beside source.
pysam numeric fetch coordinates are zero-based, half-open; FAI stores byte offsets.

```python
import json
from pathlib import Path
import pysam
p = json.loads(Path("provenance.json").read_text())
source = Path(p["source_relative_path"])
with pysam.FastaFile(str(source), filepath_index="sequences.fai") as fasta:
    records = bases = gc = long_records = 0
    for name in fasta.references:
        sequence = fasta.fetch(name)
        records += 1
        bases += len(sequence)
        gc += sequence.upper().count("G") + sequence.upper().count("C")
        long_records += len(sequence) >= 4
    assert records == p["sequence_count"] and bases == p["total_bases"]
    print(dict(records=records, bases=bases, length_at_least_4=long_records,
               literal_GC_fraction=gc / bases if bases else None))
```

Recipe identity is SHA-256 over ASCII `biov-prepared-fasta-recipe-v1` + one NUL
byte + compact UTF-8 JSON of the typed recipe in provenance field order (no
insignificant whitespace, non-ASCII unescaped). Prefix the lowercase digest with
sha256-. Output hashes are separate and are not recipe inputs. Algorithm revision
must change for any output-affecting change; schema, software package, library and
biological versions are separate concepts. There are no downloaded or inferred
provider release facts in this record.
"#,
        reference = p.recipe.reference,
        snapshot = p.recipe.snapshot_id,
        recipe = p.recipe_id,
        version = p.recipe.implementation_version
    )
}

#[cfg(test)]
mod tests;
