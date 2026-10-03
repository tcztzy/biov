BioV
====

BioV builds high-performance bioinformatics infrastructure for AI agents. Its
Python APIs provide genomic interval and sequence operations and access to
identifier-backed data; its MCP server exposes identifier discovery, provider
records and [local managed Python analysis](guides/analysis.md) with saved
results; its CLI also runs scripts directly or submits them to LSF.

BioV's capabilities work independently of an agent. Agent and workflow projects
handle task-specific tool orchestration, decisions about the next analysis step,
and biological interpretation.

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
