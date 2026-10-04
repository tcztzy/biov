# Rust datasets through MCP

This first Rust-native analysis slice keeps complete tables in Rust Polars and
uses the official Rust MCP SDK (`rmcp`) for stdio. It requires no Python runtime.
It complements the existing Python `biov mcp`; that entry point is not silently
replaced, and provider retrieval and tool-environment lifecycle are not migrated
by this slice.

The paired-export reuse extension described below adds `dataset_reopen` as the
seventh tool. Its implementation and two-process acceptance are validated in the
Linux source-built scope under SPEC T66–T67. Desktop-client integration and
portable prebuilt releases remain unclaimed.

The bounded [prepared FASTA metric tool](prepared-fasta.md#native-sequence-metric-tables) also creates native typed window datasets for these query/export/reopen tools. It requires an exact existing native snapshot and preparation, and uses strict record version 3 for its structured sequence origin; CSV exports remain version 2 and earlier records stay readable.

All BioV cache/data designs must remain quickly understandable and usable for
analysis without BioV (SPEC D0). New native exports therefore include a small
[portable result bundle](#use-a-result-without-biov), while original provider
packages may already supply suitable native documentation. This bounded slice
does not establish compliance of every existing cache; the [migration
audit](rust-migration.md#data-portability-audit-and-migration-gaps) records the gaps.

## Fresh-machine installation

Install the official Rust toolchain, using an existing Rust installation or
[rustup](https://rustup.rs/). A C toolchain is needed to build native dependencies.
The repository pins Rust 1.89.0, Polars 0.51.0 and rmcp 0.6.0; Cargo.lock fixes the
resolved dependency graph. From a checkout of this revision:

```sh
cargo test --workspace --locked
cargo install --locked --path crates/biov-cli
mkdir -p data results
biov mcp-native --data-root ./data --output-root ./results
```

Both roots must already exist. Configure the MCP client to launch that executable
with those arguments and absolute roots. Standard output is exclusively MCP
JSON-RPC; errors and lifecycle diagnostics go to standard error. Closing stdin
ends the server. Each launch creates a new in-memory dataset session.

The roots are operator-owned capability boundaries, not a filesystem sandbox
against hostile local processes. Do not let untrusted processes concurrently
replace roots, symlinks or exported files. `dataset_open` takes root-relative
`.csv` paths; `dataset_reopen` takes the relative JSON record path of a paired
BioV IPC export. Absolute paths, parent traversal and resolved symlink escapes
are rejected. Input bytes are bounded and copied once; parsing and the recorded
SHA-256 consume the same in-memory snapshot. Concurrent writes during that copy
are not an atomic source-filesystem snapshot. Exported artifacts are checked for
unchanged size and SHA-256 on retrieval.

## Complete-data workflow

For a deterministic teaching fixture, save this synthetic table as `data/counts.csv`:

```csv
sample,count
sample_a,1
sample_b,30
sample_c,20
sample_d,30
```

1. Call `dataset_open` with `{"path":"counts.csv","schema":{"count":"int64"},"preview_rows":1}`. The result
   has schema, row count 4, one preview row, unknown scientific metadata and an
   opaque session `dataset_id`. An opaque handle is not a biological accession.
2. Call `dataset_query` with that ID, `filter` equal to
   `{"column":"count","op":"ge","value":20}`, `sort` equal to
   `{"column":"count","descending":true}`, and `select` equal to
   `["sample","count"]`. The derived dataset contains all three matching rows,
   including rows absent from the opening preview. Equal sort keys retain source
   order; null sort keys go last. Execution order is filter, sort, projection.
3. Call `dataset_export` with the derived ID. The complete result is saved as an
   Arrow IPC file, strict JSON provenance record, portable manifest and README
   with unique names, without overwriting any existing output. Export returns
   its SHA-256, full row count, schema, record, execution-host paths
   (`execution_host_path`, `record_path`, `manifest_path`, `readme_path`) and a session artifact handle.
4. If the client does not share the execution filesystem, use
   `dataset_read_artifact` with that artifact ID, `offset: 0` and `max_bytes: 49152`.
   Decode each base64 chunk, append in offset order, repeat from `next_offset`
   until `eof`, and verify the complete file's SHA-256. An execution-host path
   alone is not a cross-machine transfer. This tool retrieves IPC bytes; the
   bounded JSON record is already returned by export. The README and manifest
   are saved host files, not embedded response bodies and not retrievable through
   this Arrow-only tool. Use an authorized ordinary file transfer for those
   companions when the client lacks filesystem access. Returning their host paths
   alone does not deliver a complete portable bundle to the client.
5. `dataset_release` frees a session dataset. Exported files remain on disk.

The IPC result can be read directly by Rust Polars or by Python Polars/PyArrow.
A separate Python analysis API and custom Arrow FFI are not required for that
interoperability. Python interoperability here means the standard file format,
not a promise that the current Python package wraps these dataset methods.

## Use a result without BioV

Keep the four same-directory files for one export together, preserving names:

```text
artifact_<id>.arrow          complete standard Arrow IPC file
artifact_<id>.json           mandatory strict BioV export record, version 2 for CSV / 3 for sequence metrics
artifact_<id>.manifest.json  portable descriptive manifest, version 1
artifact_<id>.README.md      ordinary-reader instructions for this result
```

Here `<id>` is the generated 32-character lowercase hexadecimal artifact ID.
Copy or move these files as a group; the original CSV or FASTA, BioV executable/package,
MCP session, repository and original absolute directory are not needed for
independent reading. The Arrow format is ordinary IPC, not a proprietary BioV
container. Existing exports with only the record and Arrow remain readable;
they do not retroactively gain companion descriptions.

Start with the generated README. It contains a directly runnable Python example
using only the standard library and `pyarrow`: load the adjacent manifest, verify
Arrow and strict-record byte counts/SHA-256, read the complete IPC table, check
its schema and full row count, then filter non-null rows and summarize numeric
columns. It does not import BioV or execute recorded lineage. New exports use the
upstream Polars writer's oldest compatibility setting: Arrow `LargeUtf8`, exposed
as `large_string` in PyArrow. The example filters the complete table directly,
without a BioV wrapper or reader-side type conversion. Other standard Arrow
readers, including Polars, can use the same file. No zero-copy claim is made. These examples demonstrate
data handling; unknown biological meanings still need scientific context before
interpreting results.

The manifest is a discoverable map of the saved result:

- `files.arrow`, `files.record`, `files.readme`: basenames relative to the
  manifest's own directory, not historical machine paths
- `content`: complete Arrow `sha256`, `bytes`, and `row_count`; `record`: strict
  record `record_version`, `sha256` and `bytes`
- `columns`: ordered `name`, `logical_type`, `polars_dtype`, `nullable`,
  `null_count`, `description`, `units` and `coordinates`. Logical types are
  `string`, `int64`, `float64` and `boolean`. `nullable: true` describes the
  allowed field schema; `null_count` reports actual missing values. Column
  descriptions/units/coordinates remain null for CSV because it has no verified
  per-column semantic declaration. Generated D15 sequence metrics supply their
  known definitions and units, plus source-relative coordinates for start/end.
  A generic CSV name such as `count` is not a definition
- `scientific_metadata`: the dataset's existing caller declarations or nulls;
  a numeric type is not a unit and an identifier does not establish sample identity
- `biological_identifier`: parsed namespace, accession, base accession and
  namespace-specific accession version when declared. Validation is syntax-only,
  not provider verification. Entry version, sequence version and provider release
  remain null; a UniProt accession does not establish them
- `source`: the recorded original CSV path, format, byte count, SHA-256 and
  input-consistency state. `historical_path` is relative to the original data
  root, explicitly not a bundle member; `required_for_reading` is false. The
  source itself is not automatically included or reverified during an independent
  result read. D15 sequence origins instead identify the registered source-relative
  FASTA path and format, with the exact native FASTA byte identity
- `lineage`: ordered operations, originally declared schema, string-safe CSV
  policy for CSV, structured `sequence_origin` for D15 metrics and any
  reopen-verification context. Its references are historical
  evidence, not files or handles needed for reading the result
- `software`, `versions`, `version_semantics`, `trust`: software facts, explicit
  unknown reference/provider versions, and the boundary between supplied claims
  and performed checks

Version numbers have separate meanings. `manifest_version: 1` versions the
companion layout; `record_version: 2` versions the strict CSV export/reopen record,
and version 3 adds the validated D15 sequence-origin schema. Record versions 1/2
remain supported without pretending to contain sequence-origin fields.
BioV and Polars versions describe software. None supplies a biological reference,
provider, entry or sequence version. A known RefSeq accession suffix may identify
that accession's namespace-defined version without resolving the table's unknown
reference assembly/context.

The manifest and README are descriptive companions. `dataset_reopen` still
requires and checks the strict record plus Arrow, and does not validate or depend
on the companions. Checksums link exact bytes to the supplied metadata; anyone
who changes the data and matching records can make them agree. This is not
producer authentication, proof of original provenance or biological QC. Treat
all stored metadata as data, not instructions to execute, and obtain a trusted
reference separately when authenticity is required.

The bundle is self-contained for understanding and reading this complete result,
not for rerunning an analysis against unbundled historical inputs. Known source
hashes and lineage remain useful even when the source paths no longer exist.
The current 64 MiB retained-data charge and all existing query limits still apply;
this does not add an out-of-core or large-cache implementation.

## Reopen a saved export after restart

Keep the exported JSON record and its IPC file together. After the first MCP
process exits, launch a second process with the old output directory as its data
root and a new output directory:

```sh
mkdir -p reopened-results
biov mcp-native --data-root ./results --output-root ./reopened-results
```

Call `dataset_reopen` with `record_path` set to the saved `.json` record's filename
relative to this new data root and optionally `preview_rows`, for example:

```json
{"record_path":"<saved-record-basename>.json","preview_rows":1}
```

Replace the placeholder with the basename of the `record_path` returned by the
earlier export. A subdirectory-relative record path is also allowed, provided
the IPC file is beside the record and both resolve within the data root. The
record's `file` field must name a basename, never an absolute path or traversal.
It must also match the recorded artifact ID and native
`artifact_<32 lowercase hexadecimal digits>.arrow` naming convention; keep the
exported IPC filename unchanged. An in-root symlink to a record resolves the IPC
beside the canonical record, rather than beside the symlink. Both resolved files
must remain within the data root.
The IPC file alone is insufficient for `dataset_reopen`; this is not an arbitrary
Arrow import mode. Independent standard readers can read the IPC file directly.
The manifest/README companions are optional for reopening and must not be passed
as `record_path`.

Reopening validates the record version and `arrow_ipc` format, file size and
SHA-256, ordered schema against the actual supported Arrow column types, full
row count and metadata fields. It parses the same bounded in-memory IPC bytes
used to compute the digest. Only the native little-endian, uncompressed, flat
string/int64/float64/boolean IPC shapes produced by this BioV slice are supported;
dictionaries, extensions, other types and unsupported Arrow metadata are rejected
before full Polars schema/array allocation. String storage may be current
`LargeUtf8` or prior `Utf8View`; other string encodings are not a general import
promise. For `LargeUtf8`, the bounded preflight checks exact 64-bit offset-buffer
length, a zero start, nonnegative monotonic offsets within the payload, valid
UTF-8/character boundaries and an end offset covering the whole values buffer.
It uses upstream Arrow metadata definitions and keeps the same allocation/file
limits before handing the complete bytes to the standard Polars reader.

The public initial-slice record version 1 remains readable. New CSV export records
use version 2 and include `reopen_verification`, which is null for a CSV-opened
dataset. Existing native `Utf8View` pairs can be reopened by this version and
re-exported as `LargeUtf8`. Earlier development readers supporting
only `Utf8View` reject the newer `LargeUtf8` outputs even though logical
record v2 is unchanged. This is one-way binary compatibility; record schema
versions do not promise older binaries understand additional physical encodings.

Success returns a fresh session dataset handle, schema, exact full row count,
bounded preview and `reopen_verification`. Use the new handle with
`dataset_preview`, `dataset_query` and `dataset_export` as usual. The entire
table is available even when only one row was previewed. The old process's
dataset and artifact handles remain expired. A subsequent export creates a new
artifact handle for bounded retrieval in the current session.

`reopen_verification` records the relative `record_path`, `record_sha256`,
`artifact_sha256`, `artifact_bytes`, `record_version` and checks for
`artifact_bytes`, `artifact_sha256`, `schema` and `row_count`. Its status fields are:

- `input_consistency`: `parsed_same_in_memory_snapshot_as_digest`
- `original_provenance`: `recorded_claims_not_independently_verified`
- `authenticity`: `not_established`

Queries retain this context and the next JSON export persists it. Reopening a
later export replaces the context with the checks of that latest pair rather
than nesting verification records indefinitely.

These are consistency checks against the supplied record, not a trusted
signature or security attestation. Someone who changes both the IPC bytes and
the record to agree can pass them. Original source history, software information
and biological metadata remain recorded claims; the checks do not authenticate
the provider, validate biology or establish QC. Missing, malformed, unsupported
or mismatched files fail explicitly, without loading the original CSV,
recomputing the result or accepting unchecked IPC. There is no persistent
catalog, daemon database, automatic restart recovery or background job.

## Contracts and limits

- UTF-8, comma-delimited CSV with a header; primitive string, int64, float64 and
  boolean columns. Unspecified columns stay strings, preserving leading-zero IDs
  and large numeric text. Use an explicit per-column `schema` for numeric/bool
  types; unknown columns, failed conversions and non-finite floats fail. CSV
  decoding uses one established `csv-core` field/record parser for preflight and
  direct typed Polars construction; Polars does not reinterpret the source CSV.
  LF, CRLF, CR and mixed record terminators are accepted. Empty physical lines
  outside quoted fields are ignored, including leading/trailing lines. Quoted
  embedded CR/LF and blank lines are preserved exactly; doubled quotes decode to
  one quote. Unterminated quoted fields, characters after a closing quote, ragged
  rows and duplicate headers fail. Quotes in a field not beginning with a quote
  are literal text under this dialect. An optional initial UTF-8 BOM is ignored
  by parsing but remains part of the recorded source-byte digest
- Unquoted empty fields become null. Quoted empty string fields remain `""`,
  distinct from null; empty numeric/boolean fields, quoted or not, become null.
  Numeric conversion permits leading ASCII spaces/tabs (whitespace-only numeric
  fields are null), but rejects trailing whitespace. Boolean `true`/`false` are
  case-insensitive. Exact int64 values never pass through float64. Reopening
  is limited to the paired native IPC exports described above. No arbitrary IPC
  import, TSV, Parquet, general FASTA import, arbitrary code, SQL, joins or aggregate operations
  in this slice
- One typed predicate (`eq`, `gt`, `ge`, `lt`, `le`, `is_null`), optional single-key
  stable sort and column projection. Null predicate outcomes do not select rows;
  use `is_null` explicitly. Boolean predicates support `eq` and `is_null`
- Integer filters remain int64 comparisons; they are not converted through
  float64. Wrong value types and missing/duplicate projection columns fail
- CSV input at most 16 MiB, 1–64 columns with names at most 128 bytes; 16 retained
  datasets and a 64 MiB conservative retained-data charge per session (Polars
  estimates plus 17 bytes per cell for view/validity overhead). This is not a
  process-memory ceiling: parsing, exports and copies require extra memory.
  Before allocating any Polars columns, CSV preflight validates all complete
  records and conversions, charging decoded string bytes, eight bytes per numeric
  cell, packed boolean/validity bounds and 17 bytes per cell against the budget
  remaining after existing datasets. Physical blank lines do not create rows or
  evade this bound. Construction uses the same immutable bytes and decoder;
  complete materialized row counts and retained charges are checked against the
  validated plan. The bounded source snapshot and field scratch buffer are
  temporary memory outside the retained-data charge
- Reopening permits a JSON record up to 64 KiB and an IPC snapshot up to 65 MiB;
  IPC footer/message metadata is limited to 1 MiB per block, with at most 4096
  record batches and 4096 buffers per batch. Record/file paths are at most 1024
  bytes. Before full IPC decoding, the allocation preflight combines logical
  payload, 17 bytes per cell and currently retained datasets against the same
  64 MiB budget. Export checks the same record/file byte caps before publishing
  files. The portable manifest is limited to 128 KiB and its README to 16 KiB;
  their caps are also checked before publication. Predictable cap failures leave
  no completed bundle. These file/parser limits do not guarantee peak memory use
- Preview defaults to 5 rows and allows 0–50. Preview row JSON has a 24 KiB budget;
  long cells are UTF-8 safely shortened to 256 bytes with truncation reported.
  Schema and metadata remain separate. Full data is preserved.
- At most 32 derivation steps and 8 KiB provenance JSON, 64 exports per session,
  and 1–49152 bytes per retrieval chunk. Each chunk rechecks the complete bounded
  artifact digest, so smaller chunks increase I/O work. Limits fail explicitly
- Metadata fields `identifier`, `species`, `reference`, `coordinates`, `units`
  default to null. Supplied values are caller declarations, not provider-verified
  science. A supported identifier is syntax-validated offline; it does not prove
  record existence or infer species, assembly, coordinates, units or QC
- Session handles cannot be reused after restart. IPC and JSON records persist
  in the configured output root and can be used as ordinary files. Explicitly
  reopening a valid saved pair creates a new handle; it does not reload an
  artifact catalog or schedule persistent jobs after restart
- Export uses create-new names. Publishing multiple files is not a directory-wide
  atomic transaction: a crash can leave a partial group of complete files.
  A new export is not reported completed until all four files have been saved.
  Storage is never automatically cleaned up

The independent `biov-identifiers` crate covers a deliberately narrow offline
RefSeq GCF/UniProt grammar. The broader Python registry/provider behavior remains
available separately. Biological ID parsing, table operations and transport are
separate crates so future CLI, Python or other consumers can reuse the libraries.

## Validation

Earlier development builds used different CSV parsers for preflight and loading.
CR-only files could lose records, and physical blank lines could add null rows.
Regenerate affected old imports from the original CSV with the corrected reader;
reopening an existing export checks its saved bytes and record, but cannot
reconstruct rows previously lost during parsing.

Portable-bundle acceptance (SPEC T68) is distinct from BioV's paired-reopen test:
move the four files away from the original data/output directories, make the
original source unavailable, and run the README's standard-reader example in an
environment where BioV cannot be imported or invoked. Verify both recorded file
hashes, ordered schema, full rows/types/nulls/duplicates/order, a filter and a
numeric summary including rows absent from previews. Check that manifest paths
remain relative, known metadata is retained, unknown semantic/version fields
stay null and no historical path is required. This gate does not certify the
unmigrated Python/provider caches. The opt-in test is
`tests/test_native_bundle_portability.py`; configure `BIOV_TEST_BINARY` with an
installed native binary and `BIOV_STANDALONE_PYTHON` with a separate environment
containing only `pyarrow==25.0.1`. The independent reader receives only four
moved files, checks that BioV is absent, and blocks attempted network/process/
SQLite access through a Python audit hook. It is distinct from the pytest driver
and does not use the MCP response as its data dictionary. The installed Linux
acceptance passed with direct filtering of the exported `LargeUtf8` table,
without casts or BioV helpers. All 50 `biov-data` tests and 17 actual MCP
subprocess cases also pass; other platforms and native client integrations
remain outside this validation claim.

`cargo test --workspace --locked` includes library contracts and real MCP stdio
subprocess tests: initialization, advertised tools, bounded open/preview, derived
complete data, IPC export and byte reconstruction/readback, invalid requests,
unknown IDs, path limits and clean shutdown. The fixture is synthetic and proves
infrastructure behavior, not a scientific analysis result.

Paired-reopen acceptance uses two real MCP subprocesses: process A
opens typed CSV, filters/sorts/selects, exports and exits on EOF; process B uses
A's output as its data root, reopens the record, previews, queries and re-exports.
Independent full IPC readback and hashes verify every expected row, type,
null, duplicate and order, including rows omitted from previews. Missing files,
tampering, schema/row-count mismatches, unsupported metadata and path escapes
return tool errors without fallback. Library cases additionally cover repeated
reopening, version 1 compatibility, malformed native IPC metadata, in-root record
aliases and allocation charging against already retained datasets. Both prior
`Utf8View` and current `LargeUtf8` cases cover empty/all-null strings, long
Unicode, multiple record batches and exact retained-budget boundaries. Malformed
`LargeUtf8` offsets/payload fail before full decoding even when the JSON digest
is updated to match the corrupted bytes.

For an installed binary outside the source directory:

```sh
BIOV_TEST_BINARY="$(command -v biov)" cargo test -p biov-cli --test mcp_stdio --locked
```

### Verified implementation scope (2026-10-03)

After the CSV parser consistency fix, the current suite passes 50 `biov-data`
tests and 130 workspace tests/doctests. All 17 real MCP subprocess cases pass against an
independently installed Linux x86_64 binary launched outside the checkout with
an empty PATH. This includes restarting between export and reuse, full-byte
reconstruction and typed readback. Adversarial metadata tests refresh the JSON
hash after mutating IPC metadata, so a matching digest cannot bypass the native
shape and allocation checks.

For the initial six-tool slice, on a fresh Linux x86_64 cloud checkout, Rust
1.89.0 and the locked dependencies were installed from official sources. The Rust
workspace tests include the
real-server MCP client flow; at that initial milestone, seven subprocess cases passed against
an independently installed binary launched outside the checkout with an empty
PATH. The complete source distribution retains all workspace members and its
lockfile; source-distribution tests, wheel build and installed sequence bindings
are separately checked. This validates a source-built Linux slice, not a portable
prebuilt release, macOS runtime or a particular desktop MCP client's integration.

The setuptools source manifest preserves the complete Cargo workspace and
lockfile, including native-only members. Creating an sdist or building a wheel
from its unpacked sources does not require Git or a checkout. The exact source
manifest and mixed native executable/extension wheel are checked separately.
