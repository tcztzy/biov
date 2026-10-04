# Rust core, MCP and optional Python

Native sequence slices are implemented: normalization, IUPAC reverse complement,
validated sequence lengths and weighted GC fractions run in a shared Rust core
with thin PyO3 batch bindings. See the [sequence
contract](sequence-contract.md) for types, errors and independent scientific
fixtures. A separate [Rust dataset MCP slice](rust-datasets.md) now opens local
CSV with string-safe schemas, executes full-data Polars queries and exports Arrow
IPC through official rmcp stdio, without Python. The implemented bounded extension,
`dataset_reopen`, reuses a saved JSON record and its paired IPC file after a
server restart; its two-process acceptance is validated in the source-built
Linux scope (SPEC T66–T67).
A separate [native-storage slice](native-storage.md) adds immutable registration
of existing RefSeq packages and explicitly declared PDB representations, with
offline filesystem discovery. The [prepared RefSeq FASTA](prepared-fasta.md) extension adds conventional FAI/TSV indices (D14/T72); its bounded `dataset_fasta_windows` tool connects exact prepared sequence access to Rust Polars tables (D15/T73). Canonical GC excludes ambiguity and is null without canonical bases; weighted GC retains the existing core IUPAC policy with every base in its denominator. Other providers and generalized prepared transforms remain separate future work.
A [native locked-tool bridge](environments.md#native-rust-locked-tool-migration)
now delegates bundled Samtools/GOATOOLS setup and literal local execution to pinned
Pixi, with recorded setup inspection and cross-task reuse. It does not complete
installed inventory, upgrades, removal or cache cleanup.
The rest of the rewrite remains planned; pandas, Biopython and RuRanges'
Rust-backed interval kernels still support unmigrated features.

Rust is a deliberate choice for AI-assisted engineering: compiler-checked types,
ownership and explicit errors can catch classes of mistakes early. This does not
establish scientific correctness, productivity gains or speedups. The tool-manager
product goal draws on `uv tool`/`uvx`; the language choice has its own rationale.

## Non-negotiable independent-data principle

Every cache and data design must let an agent quickly understand and analyze its
saved data **without BioV**. This includes provider downloads, analytical results,
existing caches being migrated and all future large-data designs. An agent must
be able to move the relevant files, identify what they contain and use an
ordinary reader without installing BioV, recovering a session or consulting the
original checkout. A standard file extension alone does not explain its meaning.

Preserve original native files. Supply the missing readable entry point, relative
file inventory, schema/data dictionary, known context, content identities and
lineage beside them; reuse adequate native provider metadata rather than replacing
it. Unknown meanings, units, coordinates, samples or versions must remain unknown.
Historical paths are provenance, not requirements for reading a moved bundle.
Manifest/record versions and software versions never stand in for biological
versions. Checksum agreement is consistency evidence, not authenticity or QC.

The current [native result bundle](rust-datasets.md#use-a-result-without-biov)
implements this for bounded Arrow exports. It does not migrate every BioV cache,
raise the native session's 64 MiB retained-data charge, or establish large-data
support. SPEC D0/V77 apply to every new design; the existing gaps below need
separate migrations and acceptance tests.

## What changes, and what must remain trustworthy

BioV is actively developed. Breaking changes to names, signatures, dataframe
classes and sequence objects are acceptable when they improve the design. Keep
Python convenient to use; do not reproduce pandas inheritance or Biopython types
just to keep old code unchanged. Document each changed contract and update its
examples in the same migration slice. No permanent dual backend is planned.

Scientific meaning and data integrity are a separate concern: preserve explicit
reference versions, coordinate conventions, units, duplicate records, meaningful
metadata and complete outputs. Revisions to scientific behavior require their own
evidence and explanation. Existing behavior is a differential reference, not a
reason to preserve a bug or an awkward API.

## Current surface to inventory before replacing it

These are actual source entry points, not promises to preserve every type:

- `biov.__init__`: exports `BioDataFrame`, `Seq`, `SequenceArray`,
  `SequenceDtype`, `SequenceValidationError`, identifier/artifact APIs,
  `read_fasta`, `read_gff3`, `align_paired_reads` and `settings`
- `dataframe.py`, `ranges.py`, `io/gff.py`: pandas subclass, overlap/intersection/
  subtraction/nearest operations, GFF import/export, arbitrary user columns and
  pandas indexes. GFF native coordinates and interval coordinates are distinct
- `seq.py`: pandas `biov.dna`, `biov.rna`, `biov.protein` dtypes and `.seq`
  accessor; the separate scalar `Seq` inherits Biopython's `Seq`
- `io/fastx.py`: `read_fasta` returns an ID-keyed dictionary of `Bio.SeqRecord`
  objects, including descriptions; record contents must not be lost merely
  because the replacement representation differs
- `artifacts.py`, `identifiers.py`, `filesystem.py`: identifier parsing,
  `Artifact`/PathLike, open handles and installed read-only fsspec schemes
- `environments.py`, `software.py`, `execution.py`, `remote.py`, `analysis.py`:
  manager setup, native argv, local/SSH/LSF boundaries, managed-analysis request
  models, persistent records, complete outputs and bounded previews
- `cli.py`: `--config`, `mcp`, `analyze`, `inspect-analysis`, `setup`, `pixi`,
  `exec`, `run`, `update`; `update` currently refreshes the identifier registry,
  not installed software
- `mcp.py`: `parse_identifiers`, `resolve_identifiers`, `run_analysis`,
  `inspect_analysis`, provider/identifier resources and registered output reads
- `crisprprimer`, `biov.azimuth`, `alignment.py`: existing Python entry points,
  score assets, native tool dispatch and BAM sort/index behavior; moving the
  dispatcher does not authorize replacing the scientific tools or models

## Data portability audit and migration gaps

This is a source audit of the paths below on 2026-10-03, not a claim that every
cache, scientific asset or output in the distribution has been exhaustively
validated. These paths already retain useful ordinary files, but do not yet all
satisfy SPEC D0. The new Rust result companions do not retroactively change them.

- **RefSeq genome packages** (`src/biov/artifacts.py`, `_refseq_gcf_path`): retain
  the complete extracted NCBI Datasets package, without flattening or renaming.
  An actual `datasets 18.38.0` download of `GCF_000005845.2`, inspected on
  2026-10-03, contains native `README.md`, `md5sum.txt`, a relative file catalog,
  `assembly_data_report.jsonl` and `sequence_report.jsonl`. These already supply
  discoverable file roles, checksums, assembly/annotation facts and sequence-ID
  mappings; do not duplicate them with a BioV replacement manifest. Moving this
  package and reading it with a Python standard-library script, with no BioV
  imports or network calls, reproduced the same results: all seven native MD5 entries and
  catalog file sizes matched; the genome contained 4,641,652 bases, with 4,300
  protein FASTA records and 4,318 CDS FASTA records. These are one real-package
  inspection's results, separate from synthetic native-table acceptance. The native
  README links [NCBI package/JSONL guidance](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/tutorials/working-with-jsonl-data-reports/)
  but has no local analysis example. A small local
  usage entry point/download receipt may fill that gap while preserving every
  upstream member. This one observed package is not proof about every provider
  package. Requested RNA was absent from the files/catalog and must be reported
  unavailable. Its assembly release (2013-09-26) and annotation release
  (2026-09-02) differ: assembly accession version does not freeze annotation
- **UniProt and other provider files** (`artifacts.py`, `_uniprot_path`,
  `_file_path`; `file_sources.py`; `ncbi_files.py`): preserve independently cached
  UniProt entry JSON and FASTA, plus native downloaded representations for other
  kinds. Format/accession validation alone does not establish retained-source
  identity or every relevant version. An AlphaFold model-selection response, for
  example, is used during retrieval but not saved alongside the chosen model.
  Inspect each provider's native metadata first and retain only missing context
  needed for independent use, without rewriting provider originals or inventing
  entry/sequence/provider versions
- **ENCODE files** (`artifacts.py`, `_encode_path`, `_encode_artifact`): preserve
  the original downloaded filename/compression and raw `metadata.json`.
  Download checks the provider's published MD5; current reuse checks metadata,
  file identity and size, without recomputing that MD5. This useful native
  evidence should be reused. Independent checksum/reuse behavior and a local
  usage path must be explicit in its migration; extra files are needed only for
  facts or guidance missing from the native material
- **HTTP filecache** (`src/biov/io/_preprocess.py`, `config.py`, and
  `src/crisprprimer/id_converter.py`): delegate cached bytes and their mapping to
  fsspec, with `BIOV_HOME` as a default unless fsspec is configured otherwise.
  They do not provide a documented independent inventory/usage path or
  a moved-directory analysis example. Document or export the source mapping and
  relevant native format/version facts; do not require BioV to recover which
  resource a cached file represents
- **Managed Python analysis** (`src/biov/analysis.py`): retains input snapshots,
  complete outputs, code, parameters, logs, environment details and SHA-256
  identities in normal files. Its `record.json`, `inputs.json`, output/snapshot
  URIs and working directory contain original absolute locations. There is no
  portable relative inventory/README contract. CSV dtypes are inferred for a
  bounded preview, not declared as a full-table schema or semantic data
  dictionary. Add a movable inventory and independent-reader acceptance without
  weakening source stability, output identity, interruption or retention checks
- **CRISPR/BLAT caches and results** (`src/crisprprimer/__init__.py`, `_blat` and
  the spacer/pair Parquet writers): use standard Parquet; BLAT caches are
  partitioned by `5-mer` and reused by query sequence within a caller-selected
  reference cache. The persisted cache does not record reference content identity,
  exact reference version, command/tool provenance or a portable dictionary for
  PSL-derived columns. Choosing a separate directory is currently the caller's
  responsibility; the cache itself cannot establish that an identically named
  query was aligned to the same reference. Add these facts and validation before
  treating such a cache as independently interpretable or safely transferable

Software/model caches and other output writers need the same per-slice inventory
when migrated; this audit is not evidence of their compliance. Native manager
locks or source checkouts can supply useful version facts but do not replace a
data-specific meaning/usage contract. Do not serialize a BioV/Python object as the
only usable data representation.

The new native-store library is an explicit registration path, not an automatic
migration of any cache audited above. Its RefSeq checks and declared PDB contract
do not make existing UniProt, AlphaFold, GEO or other legacy caches immutable or
portability-compliant. T69 tracks these open migrations. For each slice, preserve native scientific
semantics and known lineage, make unavailable metadata explicit, move a complete
bundle to a new directory, and use a standard reader with BioV absent to verify
identities/schema and analyze complete data. Test meaningful filtering or
summary, not just file existence or the first preview rows. A source checksum
records an identity; do not require an unavailable historical source or report it
as reverified merely because the result's checksum matches.

### Reproduce the real RefSeq package inspection

This opt-in download is separate from default offline tests and the native Rust
synthetic-table acceptance. The observed source is the [NCBI assembly record for
GCF_000005845.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/),
*Escherichia coli* K-12 MG1655. Use the [official NCBI CLI installation
instructions](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/).
The inspected Linux AMD64 executable came from NCBI's [official rolling binary
URL](https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets)
and reported `datasets version: 18.38.0`; that URL is not a version-specific pin.
In an empty working directory, the exact download invocation was:

```sh
datasets version
datasets download genome accession GCF_000005845.2 \
  --include gff3,rna,cds,protein,genome,seq-report \
  --filename ncbi_dataset.zip --no-progressbar
python -I -m zipfile -t ncbi_dataset.zip
python -I -m zipfile -e ncbi_dataset.zip package
(cd package && md5sum -c md5sum.txt)
```

On 2026-10-03, the ZIP was 4,154,411 bytes with SHA-256
`ca7b17b300aca5c098b3140eb648b3845a02c63d1996072e46ee0ae09684577e`.
This identifies that acquisition, not every later download: annotation or archive
metadata can change without changing assembly accession `.2`. Its native tree was:

```text
package/
  README.md
  md5sum.txt
  ncbi_dataset/data/
    dataset_catalog.json
    assembly_data_report.jsonl
    GCF_000005845.2/
      GCF_000005845.2_ASM584v2_genomic.fna
      genomic.gff
      protein.faa
      cds_from_genomic.fna
      sequence_report.jsonl
```

The standalone standard-library example is `scripts/inspect_refseq_example.py`.
Copy that one script to the download directory; it does not require the rest of
the repository or an installed BioV package. After extraction:

```sh
python -I inspect_refseq_example.py package > original.json
cp -R package moved-package
python -I inspect_refseq_example.py moved-package > moved.json
cmp original.json moved.json
```

It is deliberately scoped to this small `GCF_000005845.2` example. It is not a
general validator for other assemblies, circular features or spliced CDS.

Catalog `filePath` values resolve relative to `ncbi_dataset/data`, not `package`.
Use each `fileType` to select its representation. The report-only catalog group
has no accession; select the accession-bearing group explicitly. Preserve native
camelCase report keys: CLI summary JSON is a different representation.

The recorded inspection also verified all GFF sequence IDs/coordinates against
the genome and exact agreement of the GFF `protein_id` set with protein FASTA
IDs. The simple plus-strand `thrL` CDS at `NC_000913.3:190..255` gives 66 bases
using Python slice `[189:255]`, equal to its CDS FASTA sequence. This is a narrow
example of GFF's 1-based inclusive conversion, not a general spliced/circular-CDS
validator. Gene/pseudogene rows were 4,506/145; the 4,340 CDS feature rows are not
the number of distinct proteins. Missing RNA FASTA does not imply absent RNA
genes. A copied package produced identical analysis JSON with isolated Python
(`-I`), no BioV imports and no network calls; network access was not separately
blocked by a firewall or network namespace.

Keep these native files intact. A small acquisition receipt or local usage
example can add command/time/hash facts and explain entry points; neither needs
to duplicate the provider's biological reports or change the RefSeq cache design.
No downloaded ZIP/genome/report payload is included in this repository.

## Implementation shape

Use responsibility boundaries for offline identifiers, formats, biological
computation, provider/data/provenance, cache mechanics, tool lifecycle/execution,
optional Python and CLI/MCP. There is no fixed small-crate budget and no reason
to create empty crates. The current workspace adds `biov-identifiers`,
`biov-data`, `biov-cli`, `biov-tools` and the separate `biov-storage` library to the existing
`biov-core`/`biov-python` sequence slice.
Transport adapters call standalone Rust library APIs.

Rust-side Polars is the primary analytical engine through MCP. The dataset API
preserves identifier text by default and requires explicit numeric schemas.
Arrow IPC files provide complete typed-table interoperability with Python
Polars/PyArrow and other Arrow consumers, without a parallel Python analysis API
or custom cross-process FFI. New native exports use the upstream writer's
compatible `LargeUtf8` string encoding for direct standard-reader filtering.
This changes physical Arrow storage, not the string-safe logical schema or strict
record v2. The existing small PyO3 sequence interface remains useful without
becoming a mandatory facade for new analysis features.

Saved dataset reuse is explicit: reopening requires a BioV export record and its
same-directory IPC file under the configured data root, validates their byte
identity and table structure, and creates a new session handle. It is not
arbitrary IPC import, a persistent dataset catalog or a background-job facility.
The original handles still expire at shutdown. Reopening checks consistency
against the supplied record; a changed table and correspondingly changed record
cannot be distinguished from the original by those checks. Recorded provenance
and biological metadata are not authenticated, and reopening verification is
retained in subsequent exports. New exports also have a same-directory README
and manifest with relative file names for direct use without BioV. They do not
change the strict CSV record v2 schema or become prerequisites for paired reopening.
D15 sequence-origin exports use strict record v3 to describe the additional typed
lineage; prior CSV records v1/v2 remain readable.
The current reader accepts prior `Utf8View` and new `LargeUtf8` native strings;
old results re-export through the standard compatible writer. Earlier development readers supporting
only `Utf8View` cannot reopen new
`LargeUtf8` outputs. The bounded metadata/allocation checks were extended for
this precise physical encoding; this is not a general Arrow-import rewrite or
an increase to the 64 MiB retained-data budget. See the [dataset contract](rust-datasets.md)
for exact scope, limits and the separate acceptance status.

The new `biov-rs mcp --data-root DIR --output-root DIR` is independently tested
through the official Rust MCP SDK. The existing Python `biov mcp` retains its
separate identifier/provider/managed-analysis surface; the native slice does not
claim parity or silently replace it. Python analysis scripts and external
scientific programs continue to run in their selected environments. See the
[native dataset guide](rust-datasets.md) for current limits and installation.

## Library decisions to validate

Polars 0.51.0 and rmcp 0.6.0 are pinned and tested for the local-table slice.
Other entries below remain candidates or retained binding infrastructure; none
is a performance claim:

- [PyO3](https://pyo3.rs/) and [maturin](https://www.maturin.rs/): native Python
  extension and mixed-package distribution; test binding lifetimes, exception
  translation and cancellation, and release the interpreter only for safe
  Rust-only work
- [Official Rust MCP SDK](https://github.com/modelcontextprotocol/rust-sdk):
  used by the native dataset stdio adapter; validate tool/resource schemas,
  negotiation, errors and client behavior rather than implement a custom protocol
- [Polars](https://docs.pola.rs/): columnar computation, with explicit selection of
  features and evaluation of compile time, wheel size, memory and conversion cost
- [ruranges-core](https://github.com/pyranges/ruranges-core): evaluate the direct
  Rust kernel first, since BioV already uses RuRanges. Confirm coordinate types,
  dependency versions and operation semantics; Python-wrapper behavior is not
  automatically core behavior
- [Rust-Bio](https://docs.rs/bio/latest/bio/): useful sequence algorithms and
  [interval trees](https://docs.rs/bio/latest/bio/data_structures/interval_tree/index.html).
  An overlap tree alone does not supply BioV's entire interval algebra
- [noodles](https://docs.rs/noodles/latest/noodles/): format-specific readers and
  writers. Its [Position](https://docs.rs/noodles-core/latest/noodles_core/position/struct.Position.html)
  is 1-based; conversion to BioV's half-open intervals must be explicit

Do not treat Rust-Bio GC or an arbitrary translator/pI function as interchangeable
with BioV's current definitions. Verify the selected crate release's actual source
before pinning; online `latest` documentation is not a reproducible reference.
Remove Biopython from BioV once migrated functions no longer need it. Its age is
not evidence against its algorithms, and caller scientific environments may still
use it independently.

## Scientific and operational gates

For each slice, retain complete reference outputs and add independent evidence:

- Intervals: 0-based half-open coordinates, chromosome/strand grouping, touching
  ends, invalid/empty/overflow inputs, duplicates, stable pair/fragment order,
  nearest ties and gap-plus-one distance. Use small brute-force oracle cases
- Sequences: declared alphabets, case, null versus empty, IUPAC reverse complement,
  weighted ambiguous GC (including empty input), complete codons, every supported
  table's 64 codons, table IDs/names, ambiguous and dual-coding stops and `to_stop`
- Protein properties: canonical-residue policy, average mass in Da, water/terminus
  correction, pKa/terminal adjustments and numerical tolerance for pI; percentages
  and residue order for composition. Record constants and their provenance
- Formats and providers: descriptions/attributes, compressed and malformed files,
  exact versioned accession selection, native file bytes and coordinate conversion
- Independent data: relative inventory, ordinary-reader example, schema/dictionary,
  known and unknown semantics/versions, source identities and lineage; moved-bundle
  checksums, full-row readback and filter/summary with BioV absent (D0)
- Execution/results: literal argv and failure propagation, no cross-host local
  paths, atomic caches/records, unknown interrupted status, no resubmission,
  complete-output reuse beyond previews and detection of missing/changed files
- Interfaces: selected Python return types/exceptions, CLI help/exit codes,
  MCP schemas/error semantics and clean stdio; test intentional breaking changes
  explicitly instead of silently weakening old assertions

A self-written algorithm needs authoritative definitions/tables, hand-calculated
or independently implemented examples, differential and property tests, and
review. Tests generated from the same new implementation are not an independent
oracle. Preserve attribution/licenses for reused constants, data and source.

## Stages and distribution

1. Inventory the surfaces above and fix the new dataframe/sequence API, explicit
   breaking changes and scientific fixtures (SPEC T58)
2. Establish the core/binding/binary and one sequence slice: normalization and
   reverse complement. Validate wheels/sdist before replacing the build backend
   (T59)
3. Move remaining sequence, interval and parser logic, selecting libraries only
   after the relevant gates pass (T60)
4. Move identifiers/cache, lifecycle/execution and managed results. Implement the
   bounded tool-lifecycle additions in this core rather than duplicate them in
   Python (T56–T57, T61)
5. Finish native CLI/MCP, remove superseded Python implementations/dependencies
   and release only tested targets (T62)

Initial native release targets are Linux x86_64 and macOS ARM64, subject to actual
build/install/runtime validation. Declare a Linux libc baseline, macOS deployment
target, Rust MSRV and supported Python versions. Other platforms require their own
validation and must not be advertised prematurely. Provide a documented source
build for users outside the wheel matrix. Do not conflate BioV's package platforms
with the bundled scientific environments, which currently target Linux x86_64.

Supported prebuilt wheels require no Rust compiler. Source builds require the
specified Rust toolchain, linker and system dependencies. The sdist must contain
all Cargo/Python sources, lock/build metadata, resources, licenses, tests and docs;
build and test it after unpacking. Evaluate abi3 against real dependencies rather
than assume it covers every interpreter/ABI. Standalone CLI archives need separate
target tests and checksums. See [maturin distribution](https://www.maturin.rs/distribution.html).


## Source builds and native validation

The Cargo workspace contains eight members. The sequence pair is `biov-core`
(library and development `biov-core` binary) and `biov-python` (PyO3 extension).
The native dataset path adds `biov-identifiers` (offline biological identifiers),
`biov-data` (Polars datasets and provenance) and `biov-cli` (the `biov-rs` MCP
binary). `biov-storage` owns local immutable native snapshots and offline
resolution without Polars or transport dependencies. `biov-prepared` consumes
verified native snapshots to produce conventional FASTA indices and portable
prepared-view provenance through a pinned upstream format library. `biov-tools`
owns the bounded native locked Pixi setup/execution bridge, independent of
Polars and transport. Rust 1.89.0 is both the pinned build toolchain and declared MSRV. PyO3 is pinned to 0.26.0; the committed Cargo
lockfile controls its transitive dependencies. Maturin 1.15.0 is the pinned PEP 517 build
backend. The Python 3.12 stable ABI is selected for this small string/list-only
boundary; no NumPy ABI or interpreter objects cross detached Rust computation.
CPython 3.12–3.14 are the acceptance matrix. This does not claim free-threaded
Python or alternate-interpreter support.

For a checkout or unpacked sdist, install the official Rust toolchain, a C linker
(`cc` on Linux; Xcode Command Line Tools on macOS), Python >=3.12 and uv. Then:

```sh
cargo test --workspace --locked
uv sync --locked
uv run --locked pytest tests/ -q
uv build --wheel
```

The sdist includes both Cargo.lock and uv.lock for reproducible reference tests.
The same wheel command works from either a checkout or an unpacked sdist; install
the resulting wheel with `uv pip install <wheel>`. To create a new sdist as well
as a wheel, run `uv build` from a Git checkout with tracked sources. The configured
Git sdist generator requires that checkout, so do not use plain `uv build` to
rebuild an unpacked sdist.

`scripts/check_sdist.py dist/biov-*.tar.gz` checks exact parity with the selected
source paths in the current build checkout, as well as required files and package
boundaries. Run it against the matching, unchanged checkout immediately after
building. Adding, removing or renaming source files afterward, including untracked
files under selected source directories, changes that expected set and can report
a mismatch. This is a checkout-parity gate, not standalone inspection of an old
archive from an unrelated or later working tree.

The build needs access to crates.io and the Python package index unless their
artifacts are cached. Rebuild the extension with `uv sync --reinstall-package biov`
after Rust edits; Python-only source edits remain editable. Never add a Python
fallback to make an unbuilt checkout import successfully.

```sh
cargo run --locked -p biov-core --bin biov-core -- reverse-complement dna ACGTRYN
# NRYACGT
```

This development binary supports only `normalize` and `reverse-complement` with
an explicit kind and one sequence argument. It returns one normalized sequence
and newline on stdout, or an error on stderr and exit status 2. `--help` exits 0.
It is not bundled as a standalone release artifact and does not replace the
existing Python `biov` CLI/MCP entry points.

Local Linux x86_64 wheels are smoke-tested outside the checkout; the sdist is
unpacked, rebuilt and tested separately. Such `linux_x86_64` wheels are host-built
artifacts, not manylinux release promises. The portable Linux release gate is
manylinux_2_28 (glibc >=2.28); the macOS ARM64 gate is deployment target 11.0.
Neither target is advertised as released until a wheel built for that target
passes native install/import/runtime tests. The macOS target cannot be validated
by a Linux cross-build. CI keeps the Python acceptance matrix and Rust tests,
and verifies package assets plus wheel/source rebuild behavior. A separate
portable-wheel release workflow remains outstanding.


The local acceptance run installed the same CPython-3.12-abi3 BioV wheel under
CPython 3.12, 3.13 and 3.14. The 3.14 installation had to build the existing
`ruranges==0.2.7` dependency from source using Rust. Thus a compiler-free *complete
installation* on 3.14 is not established by BioV's abi3 wheel. The release gate
must check availability of wheels for all retained dependencies as well as BioV;
source-build tooling remains necessary when any required wheel is missing.
