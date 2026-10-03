# Rust core and Python interface

This is an accepted direction, not a description of a completed rewrite.
BioV's main logic will move to Rust and remain accessible through a thin Python
interface, a CLI and MCP. The current implementation is predominantly Python,
with pandas, Biopython and RuRanges' Rust-backed interval kernels.

Rust is a deliberate choice for AI-assisted engineering: compiler-checked types,
ownership and explicit errors can catch classes of mistakes early. This does not
establish scientific correctness, productivity gains or speedups. The tool-manager
product goal draws on `uv tool`/`uvx`; the language choice has its own rationale.

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

## Implementation shape

Use one Rust core, a PyO3 binding and a Rust binary in a small Cargo workspace.
The core owns biological operations, file/identifier handling, configuration,
cache and lifecycle state, execution and result records. Python only adapts calls,
objects, errors and necessary ecosystem hooks. Batch operations cross the binding
boundary once per batch rather than calling Python for each row.

Polars is the preferred dataframe candidate and may become the public Python
surface. Select concrete column schemas, sequence representation and null/order
rules first. Polars has no pandas-style index and distinguishes NaN from null;
conversion is not behavior-preserving by default. Unsupported object columns must
be explicitly represented or rejected. Do not silently coerce them to strings.
Benchmark complete Python-to-Rust-to-Python calls as well as the Rust operation.
See the [official pandas migration guide](https://docs.pola.rs/user-guide/migration/pandas/).

Move computation before transport. The existing Python CLI and MCP SDK can call
the Rust core during migration. Then move CLI dispatch to the Rust binary, and
move MCP to a native stdio implementation only after protocol and real-client
tests pass. A Rust binary may explicitly launch the packaged Python MCP shell
in the interim; do not call that distribution Python-free. Python analysis scripts
and external scientific programs continue to run in their selected environments.

## Library decisions to validate

These are engineering recommendations, not already selected or benchmarked pins:

- [PyO3](https://pyo3.rs/) and [maturin](https://www.maturin.rs/): native Python
  extension and mixed-package distribution; test binding lifetimes, exception
  translation and cancellation, and release the interpreter only for safe
  Rust-only work
- [Official Rust MCP SDK](https://github.com/modelcontextprotocol/rust-sdk):
  candidate for the native stdio adapter; validate tool/resource schemas,
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
