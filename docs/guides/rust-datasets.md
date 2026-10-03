# Rust datasets through MCP

This first Rust-native analysis slice keeps complete tables in Rust Polars and
uses the official Rust MCP SDK (`rmcp`) for stdio. It requires no Python runtime.
It complements the existing Python `biov mcp`; that entry point is not silently
replaced, and provider retrieval and tool-environment lifecycle are not migrated
by this slice.

## Fresh-machine installation

Install the official Rust toolchain, using an existing Rust installation or
[rustup](https://rustup.rs/). A C toolchain is needed to build native dependencies.
The repository pins Rust 1.89.0, Polars 0.51.0 and rmcp 0.6.0; Cargo.lock fixes the
resolved dependency graph. From a checkout of this revision:

```sh
cargo test --workspace --locked
cargo install --locked --path crates/biov-cli
mkdir -p data results
biov-rs mcp --data-root ./data --output-root ./results
```

Both roots must already exist. Configure the MCP client to launch that executable
with those arguments and absolute roots. Standard output is exclusively MCP
JSON-RPC; errors and lifecycle diagnostics go to standard error. Closing stdin
ends the server. Each launch creates a new in-memory dataset session.

The roots are operator-owned capability boundaries, not a filesystem sandbox
against hostile local processes. Do not let untrusted processes concurrently
replace roots, symlinks or exported files. Inputs are root-relative `.csv` files;
absolute paths, parent traversal and resolved symlink escapes are rejected.
Input bytes are bounded and copied once; parsing and the recorded SHA-256 consume
the same in-memory snapshot. Concurrent writes during that copy are not an
atomic source-filesystem snapshot. Exported artifacts are checked for unchanged
size and SHA-256 on retrieval.

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
   Arrow IPC file and a JSON provenance record with unique names, without
   overwriting any existing output. Export returns its SHA-256, full row count,
   schema, record, execution-host paths and a session artifact handle.
4. If the client does not share the execution filesystem, use
   `dataset_read_artifact` with that artifact ID, `offset: 0` and `max_bytes: 49152`.
   Decode each base64 chunk, append in offset order, repeat from `next_offset`
   until `eof`, and verify the complete file's SHA-256. An execution-host path
   alone is not a cross-machine transfer. This tool retrieves IPC bytes; the
   bounded JSON record is already returned by export.
5. `dataset_release` frees a session dataset. Exported files remain on disk.

The IPC result can be read directly by Rust Polars or by Python Polars/PyArrow.
A separate Python analysis API and custom Arrow FFI are not required for that
interoperability. Python interoperability here means the standard file format,
not a promise that the current Python package wraps these dataset methods.

## Contracts and limits

- UTF-8, comma-delimited CSV with a header; primitive string, int64, float64 and
  boolean columns. Unspecified columns stay strings, preserving leading-zero IDs
  and large numeric text. Use an explicit per-column `schema` for numeric/bool
  types; unknown columns, failed conversions and non-finite floats fail. No TSV, Parquet,
  FASTA, arbitrary code, SQL, joins or aggregate operations in this slice
- One typed predicate (`eq`, `gt`, `ge`, `lt`, `le`, `is_null`), optional single-key
  stable sort and column projection. Null predicate outcomes do not select rows;
  use `is_null` explicitly. Boolean predicates support `eq` and `is_null`
- Integer filters remain int64 comparisons; they are not converted through
  float64. Wrong value types and missing/duplicate projection columns fail
- Input at most 16 MiB, 1–64 columns with names at most 128 bytes; 16 retained
  datasets and a 64 MiB conservative retained-data charge per session (Polars
  estimates plus 17 bytes per cell for view/validity overhead). This is not a
  process-memory ceiling: parsing, exports and copies require extra memory
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
  in the configured output root and can be used as ordinary files. This slice
  does not reload an artifact catalog or schedule persistent jobs after restart
- Export uses create-new names. A crash between IPC and record persistence can
  leave an orphan complete IPC file; no completed export is reported without both
  files. Storage is never automatically cleaned up

The independent `biov-identifiers` crate covers a deliberately narrow offline
RefSeq GCF/UniProt grammar. The broader Python registry/provider behavior remains
available separately. Biological ID parsing, table operations and transport are
separate crates so future CLI, Python or other consumers can reuse the libraries.

## Validation

`cargo test --workspace --locked` includes library contracts and real MCP stdio
subprocess tests: initialization, advertised tools, bounded open/preview, derived
complete data, IPC export and byte reconstruction/readback, invalid requests,
unknown IDs, path limits and clean shutdown. The fixture is synthetic and proves
infrastructure behavior, not a scientific analysis result.

For an installed binary outside the source directory:

```sh
BIOV_TEST_BINARY="$(command -v biov-rs)" cargo test -p biov-cli --test mcp_stdio --locked
```

### Verified implementation scope (2026-10-03)

On a fresh Linux x86_64 cloud checkout, Rust 1.89.0 and the locked dependencies
were installed from official sources. The Rust workspace tests include the
real-server MCP client flow, and the same seven subprocess cases pass against
an independently installed binary launched outside the checkout with an empty
PATH. The complete source distribution retains all workspace members and its
lockfile; source-distribution tests, wheel build and installed sequence bindings
are separately checked. This validates a source-built Linux slice, not a portable
prebuilt release, macOS runtime or a particular desktop MCP client's integration.

Maturin's supported Git sdist generator preserves the full workspace instead of
pruning native-only members while retaining their lockfile entries. Creating a
new sdist from the repository requires Git and tracked sources; installing or
building a wheel from the unpacked sdist does not require a Git checkout.
