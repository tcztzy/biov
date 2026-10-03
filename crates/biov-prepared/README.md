# biov-prepared

Runtime-independent preparation of one exact, registered plain RefSeq genomic
FASTA. This crate depends on `biov-storage` and the pinned established
`noodles-fasta =0.66.0` reader/indexer/writer, not Polars, Arrow, Tokio, MCP or
Python. It does not copy, rewrite, hardlink, normalize or add files inside native
snapshots. No downloads, generic transforms, compression, GC or export command
are implemented.

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
