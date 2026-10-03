# BioV project instructions

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
  The existing Python `biov mcp` route remains separate from the in-progress
  `biov-rs mcp --data-root DIR --output-root DIR` dataset route. Do not imply
  feature parity, change the existing entry point or claim client acceptance
  until the corresponding behavior is tested
- Breaking API/type changes are acceptable during active development. Document
  them and preserve scientific semantics/data integrity rather than adding
  legacy pandas, SeqRecord or full Python-API compatibility layers
- CSV columns remain strings unless explicitly typed. Preserve leading-zero
  identifiers and numeric-looking text exactly. Preflight headers, row width,
  duplicate names and conservative cell allocation before Polars parsing;
  reject invalid declared types and non-finite Float64 values
- Operate on complete datasets, never preview rows. Keep source and derived
  dataset handles distinct, bound previews and diagnostics, confine files to
  configured roots, and independently read back complete Arrow IPC exports
  reconstructed through bounded retrieval. Record known input/output content
  identities and keep persistent exports separate from session handle lifetime.
  Session handles are not biological identifiers or durable saved-result IDs.
  Unknown reference versions, coordinates, units and sample identity stay unknown
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
planned scope separately; SPEC D12 and T63–T65 are the local-table acceptance
work until actual validation establishes support.
