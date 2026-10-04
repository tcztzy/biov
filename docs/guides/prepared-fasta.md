# Prepared RefSeq FASTA indices

This bounded slice turns a verified, pinned native RefSeq genome FASTA into a
conventional FASTA index and a readable sequence dictionary. The original native
package remains unchanged. The saved index works with ordinary indexed FASTA
readers without a running BioV process, Python BioV package or catalog database.

This is implemented and validated in the source-built Linux scope below. It does
not imply a general conversion engine, compressed FASTA support or scientific
quality control.

## Exact input and output

First [register a complete native RefSeq package](native-storage.md). Preparation
requires its canonical versioned reference, exact immutable snapshot ID and one
source-relative file selected from that snapshot's `genome_fasta` representation.
There is no implicit latest snapshot, first-file selection or remote download.

```json
{
  "reference": "refseq.gcf:GCF_000005845.2",
  "snapshot_id": "sha256-<complete source snapshot digest>",
  "source_path": "ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
}
```

Replace the placeholder with the registration result's exact snapshot ID. The
source path comes from the native catalog/registered representation; do not assume
every NCBI download uses this example's filename.

```sh
biov prepared fasta --store-root ./store --request-file prepare.json
```

The `prepared_fasta` MCP tool exposes the same request through the existing
`biov mcp-native --data-root DIR --output-root DIR --store-root DIR` server. A configured
store is required. Input parameters cannot select another store or arbitrary
external input. File contents stay on disk; the response contains bounded paths,
recipe identity, sequence counts and whether a verified preparation was reused.

```text
store/
  artifacts/refseq.gcf/GCF_000005845.2/snapshots/sha256-<source>/
    acquisition.json
    checksums.sha256
    README.md
    source/<complete native package>
  prepared/sha256-<recipe>/
    sequences.fai
    sequences.tsv
    provenance.json
    README.md
```

`sequences.fai` is the standard five-column FAI index: sequence name, length in
bases, first-base byte offset, bases per complete line, and bytes per complete
line including its terminator. `sequences.tsv` has a header and two columns:
`sequence_id` (exact opaque string, including leading zeroes) and `length`
(unsigned base count). Rows preserve source record order. It is a sequence
dictionary, not another copy of the genome. Descriptions remain in the FASTA.
The generated README gives the exact schema and ordinary-reader command. Sequence identifiers, case and
native FASTA bytes are preserved. This operation does not infer genes, taxonomic
identity, reference quality, topology or annotation correctness from an accession.
The native package's catalog and reports remain authoritative for their metadata.

## Upstream indexing and supported input

The Rust library uses the established
[noodles-fasta 0.66.0 indexer](https://docs.rs/noodles-fasta/0.66.0/noodles_fasta/io/struct.Indexer.html)
and its FAI writer. This version supports the pinned Rust 1.89 toolchain. BioV
validates the supported input envelope and bounds memory before passing complete
physical lines to the upstream parser; it does not implement an alternative
FASTA offset or wrapping algorithm.

The initial input contract is deliberately narrow: plain uncompressed genomic
FASTA, printable ASCII headers (horizontal tab is also allowed as a separator), unique nonempty identifiers, case-preserving
IUPAC DNA symbols and LF or CRLF line endings. Bare CR, embedded sequence
whitespace, blank sequence lines, unsupported symbols and malformed wrapping
are rejected rather than normalized. The identifier is the upstream FASTA name
token, not the entire description. A physical line is bounded to 1 MiB including
its terminator. This limits unusually long unwrapped records, not total genome
size. Full sequences are never accumulated for indexing.

The I/O buffer is 64 KiB. Identifiers are at most 4,096 bytes and there are at
most 100,000 records. The retained duplicate-ID charge (identifier bytes plus
128 bytes per identifier) is at most 16 MiB; each generated FAI/TSV is separately
bounded to 16 MiB. Provenance and the compact response are each at most 64 KiB. These are allocation/accounting
limits, not a measured process-wide peak-memory or speed claim. Preparation and
reuse verify native snapshot integrity and scan the selected FASTA. This first
proof does not promise constant-time cache hits for large genomes.

FAI offsets and line widths describe bytes in the original file. Sequence lengths
count bases. Standard reader APIs may use different interval conventions: pysam
`fetch(name, start, end)` uses zero-based half-open coordinates, whereas samtools
region strings normally use one-based inclusive coordinates. Neither convention
is silently substituted for the other.

## Recipe, integrity and reuse

The recipe key covers exact input identity, selected source path, the indexing
implementation contract, upstream library version and output-affecting parameters.
Output SHA-256 values separately identify the actual saved bytes. A changed tool,
algorithm contract, parameter or input must produce a different recipe identity;
algorithm changes require an explicit contract-version update.

Only fully written and verified preparations are published, using atomic
no-replace semantics. Concurrent duplicate work preserves and verifies the
winner. Partial staging directories do not count as available preparations.
Repeated requests perform a complete native snapshot verification, replay the
upstream indexer into hash-only sinks, compare the regenerated output identities
and full typed provenance, rehash saved files, then verify the native snapshot
again before returning `reused: true`; corruption is an error, never an instruction to overwrite a
snapshot or silently repair historical results. Original import directories and
the registered native snapshot are not modified.

Provenance records registration/snapshot linkage and the actual indexing recipe.
Unknown acquisition facts remain unknown. SHA-256 agreement establishes local
consistency, not producer authenticity or independent biological validation.

## Moving the dependency closure

Existing root documentation is preserved, including earlier registration-only
README text. Each prepared view's README describes its current files and dependency
closure; an older root README is not an inventory of subsequently prepared views.

Raw data is reused by relative reference rather than copied for every prepared
view. Consequently the prepared directory alone is **not self-contained**. Copy
both the exact referenced native snapshot and the prepared directory, preserving
their paths relative to a common store root. Copying the complete store also
preserves this closure. Do not replace those copies with links to the old root.

The generated README and provenance provide the exact dependency path. A reader
must resolve it relative to the prepared directory, never relative to its current
working directory or an old execution-host path. No database, original download
location, live MCP handle or network connection is required. Snapshot identity
does not pin annotation bytes merely by assembly version: the full source digest
does.

Use the generated ordinary-reader example with an explicit external index path.
For example, pysam's `FastaFile` accepts `filepath_index`; it need not create an
index beside the immutable raw FASTA. This avoids changing provider-native files.

## Native sequence metric tables

The `dataset_fasta_windows` native MCP tool reads one exact sequence from an
already existing, verified preparation. Start the same native server with
`--data-root DIR --output-root DIR --store-root DIR`. It requires all three
preparation selectors plus the exact `recipe_id`, exact opaque `sequence_id`
and a positive `window_size` in bases. `preview_rows` defaults to 5 and is at
most 50. No sequence or recipe is guessed, and a missing preparation is an error;
call `prepared_fasta` explicitly first. This tool does not download, normalize or
write source/prepared files.

```json
{
  "reference": "refseq.gcf:GCF_000005845.2",
  "snapshot_id": "sha256-<complete source snapshot digest>",
  "source_path": "ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna",
  "recipe_id": "sha256-<exact preparation recipe digest>",
  "sequence_id": "NC_000913.3",
  "window_size": 10000,
  "preview_rows": 5
}
```

Use the exact identities returned by registration and preparation, not the
placeholders. Sequence identifiers are case-sensitive first-header tokens;
leading zeroes, punctuation and descriptions are not rewritten. Windows are
non-overlapping and ordered from the beginning of the selected sequence. They
use source-sequence-relative zero-based half-open `[start,end)` coordinates.
The final row is short when length is not divisible by window size; it is retained.

The table has these ordered columns:

| Column | Type | Meaning |
| --- | --- | --- |
| `sequence_id` | string | Exact selected FASTA name token |
| `start` | int64 | Inclusive zero-based start, in bases |
| `end` | int64 | Exclusive end, in bases |
| `length` | int64 | Bases in this window, `end-start` |
| `is_full_window` | boolean | Whether `length == window_size` |
| `canonical_base_count` | int64 | Case-insensitive A/C/G/T count |
| `gc_base_count` | int64 | Case-insensitive literal G/C count |
| `gc_fraction` | float64, nullable | G/C count divided by canonical count |
| `weighted_gc_fraction` | float64 | Existing IUPAC-weighted GC divided by all bases |

Canonical GC excludes every ambiguous symbol from both numerator and denominator.
An all-ambiguous window has `gc_fraction: null`, rather than zero. Weighted GC
uses the existing Rust sequence-core policy: each IUPAC symbol represents an
equal-weight base set, so N contributes 1/2, B/V 2/3, D/H 1/3, S 1 and W 0.
Every symbol counts in the weighted denominator. Both fractions are dimensionless;
case does not change either metric. For example, `NNNN` has canonical counts zero,
canonical GC null and weighted GC 0.5.

The response is the ordinary bounded dataset summary with an opaque session
`dataset_id`, typed schema, row count and preview, plus `sequence_summary`.
That summary contains `sequence_id`, whole `length`, `canonical_base_count`,
`gc_base_count`, nullable `gc_fraction` and `weighted_gc_fraction`, accumulated
in the same pass over all selected bases. Fractions are derived from whole counts,
not averages of window percentages. This whole-sequence summary is a response;
the persisted result is the complete window table and its lineage.

Upstream noodles indexed queries retain at most one requested window of sequence
bases, capped at 1 MiB. The implementation verifies the exact native and prepared
bundle before and after reading. Table row/cell/identifier memory is conservatively
preflighted against the remaining 64 MiB session retained-data charge before
allocating vectors or Polars columns. The normal 16-dataset session limit still
applies. This is neither a process-wide peak-memory guarantee nor unbounded,
lazy/disk-backed sequence analytics. Deep integrity verification may scan complete
source files even though indexed extraction retains bounded windows.

Reuse the returned handle with `dataset_query`, `dataset_export` and
`dataset_release`. Complete exports use the same four-file Arrow IPC bundle as
[other native datasets](rust-datasets.md#use-a-result-without-biov), with standard
LargeUtf8 strings and known metric column meanings/units/coordinates. Their
`provenance.sequence_origin` records the exact reference/snapshot/recipe, selected
sequence and length, window size, FAI/dictionary SHA-256 values, coordinate and
GC policies, algorithm name and revision. Native FASTA byte size/SHA-256 remain
separate source-content identities. FASTA origins do not claim a CSV schema
policy. Their mandatory strict export record is version 3; CSV exports remain
version 2, and prior record versions 1/2 remain supported by paired reopening. Query operations append lineage without dropping sequence origin.

The four-file Arrow export is self-contained for analysis of the saved windows;
it does not need the native/prepared dependency closure. That closure remains
necessary for independently recomputing metrics from raw sequence. A new native
process can `dataset_reopen` the mandatory record/Arrow pair without the store,
restoring a fresh session handle and preserving origin/operation records. Reopen
checks saved consistency, not source availability, producer authenticity or
independent scientific provenance. Unknown species/provider releases remain
unknown.

The installed acceptance gate is `tests/test_fasta_metrics_portability.py`.
It exercises tiny LF/CRLF mixed-case/IUPAC fixtures, exact names, whole summaries,
complete query/export/new-process reopen, recipe/width/sequence errors, missing
preparations and preflight rejection without consuming handle slots. It moves
full and filtered exports, removes the original source/store/output locations,
and checks full rows, hashes, schema, filtering and sums with stock
`pyarrow==25.0.1` in a separate environment without BioV and with kernel network
calls denied. `BIOV_STANDALONE_PYTHON` selects that reader. The opt-in existing
real E. coli package additionally expects 465 windows at 10,000 bases: 464 full
windows and a final 1,652-base row, totaling 4,641,652 bases and 2,357,528 literal
G/C bases. All four installed Linux acceptance cases pass (SPEC T73). This validates the
stated input/platform scope, not a general large-genome analytics claim.

## Validation and remaining scope

Executed source-built Linux acceptance covers 25 prepared-library tests, all
21 installed CLI/MCP subprocess cases and 10 independent prepared-portability
cases. The latter include the actual downloaded RefSeq package, tiny LF/CRLF
fixtures, exact leading-zero identifiers/case, complete records and boundary/
interior subsequences, recipe separation, corruption rejection and SIGKILL/retry.
The complete dependency closure is moved, original source/store paths disappear,
and the independent reader runs without BioV or database access and with network
system calls denied in that test child. Native originals remain byte-identical.

Independent adversarial review additionally checked 31 valid fixtures, 90 full
records and 4,840 indexed slice comparisons; 35 malformed-input cases and nine
integrity/path/response-budget cases were safely rejected. Coverage includes
physical-line/identifier boundaries, quoted and Unicode source filenames, forged
output-plus-provenance hashes, symlinks, extra output members and an oversized
escaped-host-path response that must fail before publication. This is test
evidence for the stated input/OS scope, not universal FASTA/platform acceptance.

The independent reader environment uses `pysam==0.24.1`, `biopython==1.88` and
`numpy==2.5.3`, with BioV absent. Default fixtures require no provider download.
The existing real GCF_000005845.2 package is opt-in through
`BIOV_TEST_REFSEQ_PACKAGE`; the reader interpreter is selected with
`BIOV_STANDALONE_FASTA`, and `BIOV_TEST_BINARY` selects the installed native CLI.

Compressed/BGZF FASTA, new provider adapters, GEO matrices, generalized prepared
transforms, automatic downloads, garbage collection and distributed storage need
separate contracts and acceptance. A conventional index is one useful
analysis-ready representation; this does not require converting all biological
formats into Arrow.

## Format references

- [FAI file conventions](https://www.htslib.org/doc/faidx.html)
- [samtools faidx and duplicate-name limitations](https://www.htslib.org/doc/samtools-faidx.html)
- [pysam indexed FASTA API](https://pysam.readthedocs.io/en/latest/api.html#fasta-files)
- [Pinned noodles indexer source](https://github.com/zaeleus/noodles/blob/e54dcd26cb6c734e2b53cc0bfb58754240fbd13a/noodles-fasta/src/io/indexer.rs)

Samtools/HTSlib was evaluated as the alternative production mechanism and is used
through pysam for independent acceptance. The Rust library avoids requiring an
external indexing executable. Upstream samtools duplicate-name warnings are not
sufficient for this preparation contract, which rejects duplicate identifiers.
