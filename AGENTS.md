# BioV project instructions

## Project scope

BioV aims to own the lifecycle of biological tools and their environments across
tasks, following the install/run distinction of `uv tool` and `uvx`. Existing
package managers provide resolution and installation. Biological input/output
contracts, reusable data access, interval/sequence APIs and execution remain
part of BioV and work without an agent or model. Distinguish the planned lifecycle
in SPEC D6–D8 from currently implemented interfaces.

- Keep task-specific tool orchestration, analysis decisions, and biological
  interpretation in agent or workflow projects that use BioV.
- Expose BioV capabilities through Python, CLI, or MCP where each interface is
  useful. The MCP server provides identifier discovery, provider records and
  local managed Python analysis with saved results; it is not the implementation
  boundary for every Python API.
- Use established analysis libraries and tools for their scientific algorithms.
  Add BioV code for reusable data contracts, access, and execution needs.
- Treat high performance as a design goal. Describe the Rust-backed interval
  implementation accurately, and make quantitative speed claims only when
  supported by reproducible benchmarks.
- Resolve local cache locations through `BIOV_HOME`, fsspec configuration, or
  explicit directory arguments. Do not commit workstation-specific cache paths,
  including in documentation examples or packaged metadata.

Keep public behavior consistent with `SPEC.md` and the documented coordinate,
sequence, identifier, and artifact contracts.
