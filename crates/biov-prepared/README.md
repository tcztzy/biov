# biov-prepared

Runtime-independent preparation of one exact, registered plain RefSeq genomic
FASTA. This crate depends on `biov-storage` and the pinned established
`noodles-fasta =0.66.0` reader/indexer/writer and indexed queries, with
`noodles-core =0.20.0` coordinates and `biov-core` nucleotide GC computation,
not Polars, Arrow, Tokio, MCP or Python. It does not copy, rewrite, hardlink,
normalize or add files inside native snapshots. No downloads, generic transforms,
compression, garbage collection or export command are implemented.

## API

```rust,no_run
use biov_prepared::{PrepareFastaRequest, PreparedStore};
let store = PreparedStore::new("existing-store")?;
let result = store.prepare_fasta(PrepareFastaRequest {
    reference: "refseq.gcf:GCF_000005845.2".into(),
    snapshot_id: format!("sha256-{}", "0".repeat(64)), // replace with actual pin
    source_path: "ncbi_dataset/data/GCF_000005845.2/genome.fna".into(),
})?;
println!("{}", result.fai.execution_host_path);
# Ok::<(), biov_prepared::PreparedError>(())
```

The reference must be canonical and biologically versioned; the pin must be the
full source snapshot digest. The source-relative path must be in the exact
snapshot's verified `genome_fasta` result. Construction only reads an existing
trusted UTF-8 store root. Preparation errors do not replace published results.
The serialized result is at most 64 KiB and contains recipe identity, reuse flag,
sequence/base counts, and store-relative plus execution-host paths to the source,
FAI, dictionary, provenance and README. Host paths are not transfer mechanisms.

### Read-only indexed window metrics

`PreparedStore::preflight_fasta_window_metrics(FastaWindowMetricsRequest)` selects
one exact sequence from an **existing** preparation. The request has flat fields
`reference`, `snapshot_id`, `source_path`, `recipe_id`, `sequence_id` and
`window_size`; the first three are the same exact native selection as preparation.
The full lowercase `sha256-` recipe pin must match the current indexing recipe.
Sequence names are exact opaque identifiers, including case, leading zeroes and
punctuation. A text such as `0001:alt` is a name, not a region expression.

Preflight independently verifies the complete native snapshot, regenerates FAI
and dictionary identities into hash-only sinks, checks the complete saved
provenance/output bundle, decodes the verified FAI using noodles and returns a
`FastaMetricsSelection`. It never creates or repairs an absent/corrupt preparation.
Selection metadata includes the exact sequence length and `window_row_count`
before a consumer allocates any result columns, plus native/recipe/source/index
identities. Its `source_path` is store-relative, including the pinned native
snapshot path; the request's `source_path` is relative to that snapshot's source.

`selection.stream_windows(|row| ...)` emits `FastaWindowMetric` rows in sequence
order. Coordinates are zero-based, half-open: `[0,window_size)`, then contiguous
nonoverlapping windows, with a final shorter window if necessary. `length` is
`end - start`, and `is_full_window` explicitly identifies full-width windows.
The positive width is at most 1,048,576 bases. Each noodles indexed query has
explicit bounded endpoints; no unbounded full-contig query is used. The selected
FAI record is retained alone during streaming, avoiding a full-index name search
for every window. LF/CRLF and native ASCII case do not change scientific counts.

Rows and the same-pass `FastaMetricsSummary` contain:

- `canonical_base_count`: the count of literal A/C/G/T, case-insensitive
- `gc_base_count`: the count of literal G/C, case-insensitive
- `gc_fraction`: G/C divided by canonical A/C/G/T only, null when none exist
- `weighted_gc_fraction`: equal-weight IUPAC GC probability, divided by all bases

For weighted GC, G/C/S contribute 1; A/T/W contribute 0; R/Y/K/M/N contribute
1/2; B/V contribute 2/3; D/H contribute 1/3. Whole-sequence aggregation reuses
`biov-core::sequence::nucleotide_gc_counts` and adds exact sixths before final
floating-point division. It is not an average of per-window fractions and does
not read the sequence a second time for its summary. All-ambiguous windows keep
their bases and weighted values rather than disappearing or becoming canonical.
The canonical metric differs from definitions that count S or W as canonical.

The input snapshot and prepared output identities are verified again after
emission and before success. A callback consumer must discard every emitted row
if streaming returns an error; a table consumer must not publish its handle
until success. Source and prepared files stay unchanged. The trusted-local-writer
boundary and consistency, rather than authenticity, caveats below still apply.

The indexed sequence buffer is capped per query; normalization for validated
core metrics uses an additional at-most-window-sized temporary string. FAI
metadata has the existing separate 16 MiB bound, and I/O buffers are 64 KiB.
These are independent allocation bounds, not an aggregate peak-memory guarantee.
This crate does not retain the emitted table or impose a total sequence byte
cap. Table/session consumers must independently preflight the complete row count
and retained-data budget; a large FASTA does not imply an unbounded table route.

## Layout and portable dependency closure

```
store/
  artifacts/refseq.gcf/<versioned-accession>/snapshots/<exact-source-id>/
    acquisition.json
    checksums.sha256
    README.md
    source/<complete native package>
  prepared/<recipe-id>/
    sequences.fai
    sequences.tsv
    provenance.json
    README.md
```

Copy the complete exact source snapshot and the complete prepared directory,
preserving both paths below any new common root. A prepared directory alone is
incomplete. The source and full snapshot are identified relative to the prepared
directory in `provenance.json`; no historical machine path is persisted. The
complete native package and its own metadata remain authoritative. Existing store
root documentation is never rewritten by this crate.

`sequences.fai` uses the conventional five columns: name, sequence length in
bases, zero-based first-base byte offset, bases per full line, bytes per full line
including its newline. `sequences.tsv` is UTF-8 TSV with header
`sequence_id\tlength`; rows contain the exact opaque identifier and unsigned
base count in source record order. Leading zeroes, punctuation and case remain
significant. There are no feature coordinates or invented chromosome roles,
annotation versions, acquisition facts, scientific units or sample assignments.
The generated README includes a standard-reader example using
`pysam.FastaFile(source, filepath_index="sequences.fai")` explicitly. This avoids
writing an index alongside immutable native data. pysam numeric query intervals
are zero-based and half-open; FAI offsets are byte offsets.

## Recipe contract

`FastaRecipe` is a typed record containing schema version, algorithm name,
algorithm revision, implementation package/version, library name/version, exact
canonical reference, source snapshot/content identity, selected source path,
selected file SHA-256/size and explicit fixed `FastaParameters`.

The recipe ID is `sha256-` followed by lowercase SHA-256 over:

1. ASCII `biov-prepared-fasta-recipe-v1` and one NUL byte
2. Compact UTF-8 JSON of `FastaRecipe`, in the declared Rust field order visible
   in provenance, with no insignificant whitespace and non-ASCII unescaped

Output hashes are separate; they do not enter the recipe hash. The algorithm
revision MUST change for any output-affecting implementation/README/validation
contract change, even if the package version is unchanged. Library version,
algorithm revision, provenance schema, biological accession version and native
snapshot digest are distinct. The public recipe type describes this fixed recipe;
it is not a user-configurable transform engine.

`FastaProvenance` schema 1 contains that exact recipe and its key, source-relative
file/snapshot paths, counts and three output identities (`sequences.fai`,
`sequences.tsv`, `README.md`, each with byte size and SHA-256). It contains no
clock-dependent field or absolute source location. Provenance JSON is not
self-hashed; reuse reconstructs the complete expected typed record independently
from the verified source and current fixed recipe.

## Supported input and bounds

- Plain uncompressed FASTA, LF or CRLF, final newline optional
- Nonempty records and unique nonempty first-whitespace-delimited identifier
  tokens immediately after `>`; descriptions remain in source
- Printable ASCII definitions, with horizontal tab additionally permitted
- IUPAC DNA ASCII `ACGTRYSWKMBDHVN` in either case, preserved exactly
- No blanks, sequence whitespace, gaps, internal `>`, bare CR, non-ASCII or other
  symbols; unsupported data is rejected, never silently normalized
- Upstream checks FAI-compatible wrapping: nonfinal lines in each record have
  equal width/base count; the last line may be shorter

Input is fed to noodles as complete bounded physical lines. This both caps its
header allocation and removes chunk-boundary ambiguity around CR and `>`.
BioV does not implement another FASTA name decoder, offset calculation, wrapping
algorithm or FAI writer.

Bounds: 64 KiB I/O buffers; one physical line up to 1 MiB including terminator;
4,096-byte identifiers; 100,000 records; retained identifier charge of identifier
bytes plus 128 per record at most 16 MiB; each FAI and TSV at most 16 MiB;
provenance, README and public response each at most 64 KiB. These bounds are
independent, not an aggregate peak-memory guarantee. There is no total sequence
byte cap. All source hashing and sequence indexing are sequential and streaming.

## Integrity and publication

Both new preparation and reuse perform initial full `NativeStore.resolve`
verification. The selected open input stream is simultaneously hashed and indexed;
its bytes and digest must match the verified receipt. Before publishing or
returning reuse, another full native resolution detects ordinary in-flight edits.
Reuse regenerates FAI and dictionary into hash-only sinks, compares the complete
expected provenance record, rehashes all saved outputs and rejects extra members.
Changing both an output and its recorded hash cannot forge a reusable result.
Missing or corrupt input/output is an error, never silent replacement or repair.
This first implementation deliberately trades repeated sequential work for
integrity; it does not claim metadata-only reuse or large-store query performance.

Writes stage under `.staging/prepare-fasta-*`, sync files and directories, acquire
a store-local per-recipe lock, and publish with atomic no-replace rename. Existing
verified winners are reused. Response-size validation occurs before publication.
An interrupted stage is invisible and a fresh process may safely retry; stale
stages are left for explicit operator cleanup, not automatic deletion/GC.

The root is a trusted operator boundary, matching native storage. Symlinks,
special files, unsafe relative paths and unexpected prepared members are rejected,
but this is not a sandbox against hostile concurrent writers or administrator
mutation. Linux is the tested platform. Other implementation branches do not
establish filesystem/platform acceptance. Checksums establish consistency rather
than provider authenticity or independent scientific QC.

Run `cargo test -p biov-prepared --locked --offline` for core tests. The repository
acceptance suite additionally relocates the exact closure and verifies complete
records and random-access slices using isolated pysam plus Biopython with BioV
absent and network disabled, including an opt-in real RefSeq package.
