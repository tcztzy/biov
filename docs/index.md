BioV
====

BioV aims to be a biology-focused tool manager: install or run a tool, reuse its
environment across tasks, and know where it lives and how to maintain it. Its
direction follows `uv tool` and `uvx`, while retaining the biological input/output
contracts needed to use scientific tools correctly together.

Today BioV provides locked scientific environments and on-demand execution,
identifier-backed data, interval/sequence APIs and local managed analysis with
saved results. A complete tool lifecycle, including unified inventory, upgrades,
uninstall and cleanup, is [planned](guides/environments.md#planned-tool-lifecycle).
These capabilities work without an agent; task-specific orchestration, method
selection and biological interpretation belong to the caller or its skills.

## Highlights

- **Agent-facing interfaces**: MCP tools for identifier discovery, provider records and local managed analysis, with bounded previews and access to complete saved results
- **Pydantic-powered**: Built-in validation and serialization for robust data handling
- **Pandas ecosystem**: Developer-friendly DataFrame operations with extended bioinformatics capabilities
- **RuRanges interval kernel**: Stable BioDataFrame range semantics over NumPy/Rust kernels
- **Typed sequence Series**: Explicit nullable DNA, RNA, and protein dtypes with a `.seq` API
- **Identifiers.org MCP**: Identifier discovery, native RefSeq/UniProt records, file representation descriptions
- **Analysis files**: Identifier-backed sequences, structures, articles, variant records, expression matrices, and other native provider files through BioV and fsspec
- **Portable analysis**: Executor-local identifier artifacts used directly by ordinary ecosystem libraries
- **Modern tooling**: Full type hints support and configuration through environment variables

BioV requires Python 3.12 or newer. Read the [genomic range contract](guides/ranges.md), [typed sequence Series contract](guides/sequences.md), [identifiers.org MCP guide](guides/identifiers.md), and [identifier-backed analysis guide](guides/artifacts.md) before using these APIs.
