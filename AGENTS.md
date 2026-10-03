# BioV project instructions

## Non-negotiable data portability principle

Every cache and data design must let an agent quickly understand and analyze the
saved data **without BioV**. This applies to provider caches, derived results,
managed-analysis outputs, and every future storage or large-data design. BioV
may make discovery and reuse easier; it must not become the only reader or the
only source of the meaning of saved data.

- Keep complete data in documented standard/native formats. Preserve provider
  bytes and native package layouts; add descriptions beside them rather than
  rewriting originals or hiding the only usable representation in opaque state
- Make a copied or moved data bundle understandable from its files: a readable
  entry point, relative inventory, format/schema or data dictionary, known
  scientific metadata, source content identities, and transformation lineage.
  Use a provider's suitable native documentation where available, with a small
  companion only for gaps. Never require a live MCP session, BioV installation,
  repository checkout, catalog database, or historical source path to read it
- Separate known facts, caller declarations, unknown and not-applicable fields.
  Never invent column meanings, units, coordinates, sample identity, reference
  assemblies or versions. Format/manifest versions, software versions and
  namespace-specific biological versions are different concepts
- Include an ordinary-reader example and an acceptance test that moves the
  bundle, removes BioV from the reader environment, verifies byte identities and
  schema, reads complete records and performs a meaningful filter/summary. Checksums
  establish consistency with the supplied record, not authenticity or scientific QC
- Existing storage that does not yet meet this rule is an explicit migration
  gap, not an exception for new designs. Keep the audit in
  `docs/guides/rust-migration.md` truthful; a portable Rust result bundle does not
  establish compliance of every Python/provider/cache path

## Project scope

BioV aims to own the lifecycle of biological tools and their environments across
tasks, following the install/run distinction of `uv tool` and `uvx`. Existing
package managers provide resolution and installation. Biological input/output
contracts, reusable data access, interval/sequence computation and execution
remain part of BioV and work without an agent or model. Keep namespace,
accession and version semantics explicit: biological identifiers are a
first-class capability, not incidental strings. Distinguish the planned
lifecycle in SPEC D6–D8 from currently implemented interfaces.

- Keep task-specific tool orchestration, analysis decisions and biological
  interpretation in agent or workflow projects that use BioV
- Put BioV-owned main logic in Rust. Rust-side Polars exposed through MCP is the
  primary analytical interface; use Arrow IPC for complete typed-table
  interoperability with Python and other consumers. Do not build a parallel
  Python analysis API or require compatibility wrappers. The existing small
  PyO3 sequence binding may remain where useful
- Follow the logical boundaries in SPEC D9: offline identifiers; formats;
  biological computation; provider resolution, data and provenance; cache;
  manager backends, environments and execution; optional Python; CLI/MCP.
  Boundaries may be modules or crates. There is no fixed tiny crate budget;
  create a crate when real responsibilities or dependencies justify it, not to
  fill an architectural diagram
- Keep CLI/MCP adapters thin and use the official `rmcp` SDK for native MCP.
  The existing Python `biov mcp` route remains separate from the validated local
  `biov-rs mcp --data-root DIR --output-root DIR` dataset route. Do not imply
  feature parity, change the existing entry point or claim client acceptance
  until the corresponding behavior is tested
- Breaking API/type changes are acceptable during active development. Document
  them and preserve scientific semantics/data integrity rather than adding
  legacy pandas, SeqRecord or full Python-API compatibility layers
- CSV columns remain strings unless explicitly typed. Preserve leading-zero
  identifiers and numeric-looking text exactly. Preflight headers, row width,
  duplicate names, declared types and conservative payload/cell allocation
  against the remaining session budget before constructing Polars columns. Use
  the same established CSV decoder for validation and materialization; preserve
  quoted empty strings separately from missing fields and reject malformed quotes
  and non-finite Float64 values
- Operate on complete datasets, never preview rows. Keep source and derived
  dataset handles distinct, bound previews and diagnostics, confine files to
  configured roots, and independently read back complete Arrow IPC exports
  reconstructed through bounded retrieval. Record known input/output content
  identities and keep persistent exports separate from session handle lifetime.
  Session handles are not biological identifiers or durable saved-result IDs.
  Unknown reference versions, coordinates, units and sample identity stay unknown
- Use the upstream Arrow writer's compatibility setting for interoperable native
  exports: new strings use LargeUtf8 for direct standard-reader filtering. Keep
  previous native Utf8View exports readable through the bounded paired-reopen
  route; do not turn its metadata/allocation preflight into a general Arrow reader
- Reopen saved native IPC through its mandatory record, with bounded metadata
  preflight before decoding. Restore new session handles rather than reviving
  expired ones. Matching bytes/schema/digests establish consistency with the
  supplied record, not producer authenticity or independently verified provenance
- Distinguish retained-data budgets from peak-memory guarantees. Charge Polars
  estimated size plus the documented cell overhead, bound query provenance, and
  require existing UTF-8 trusted roots; canonical path checks are not a sandbox
  against hostile concurrent filesystem writers
- Prefer established scientific libraries/tools. BioV-owned Rust implementations
  for uncovered algorithms require independent scientific checks and provenance,
  not merely compilation or agreement with AI-generated expectations
- Treat high performance as a design goal. Describe implemented Rust and
  Rust-backed features accurately; make zero-copy or quantitative speed claims
  only when supported by reproducible measurements
- Resolve local cache locations through `BIOV_HOME`, fsspec configuration or
  explicit directory arguments. Do not commit workstation-specific cache paths,
  including documentation examples or packaged metadata

Keep public behavior consistent with `SPEC.md` and the documented coordinate,
sequence, identifier and artifact contracts. Report current, in-progress and
planned scope separately; SPEC D0 applies throughout. D12 and T63–T68 state the
local-table, paired-record reopen and portable-result scope, with validation and
platform/client limitations explicit.
