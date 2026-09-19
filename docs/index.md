BioV
====

Next-generation development experience for computational molecular biology.

## Highlights

- **LLM friendly/native/driven**: Designed for seamless integration with large language models, built for LLM workflows, and optimized for LLM-assisted development
- **Pydantic-powered**: Built-in validation and serialization for robust data handling
- **Pandas ecosystem**: Developer-friendly DataFrame operations with extended bioinformatics capabilities
- **RuRanges interval kernel**: Stable BioDataFrame range semantics over NumPy/Rust kernels
- **Typed sequence Series**: Explicit nullable DNA, RNA, and protein dtypes with a `.seq` API
- **Identifiers.org MCP**: RefSeq/UniProt data resources, generic registry resources, prompt parsing, and raw-response synchronization
- **Portable analysis**: Executor-local identifier artifacts used directly by ordinary ecosystem libraries
- **Modern tooling**: Full type hints support and configuration through environment variables

BioV requires Python 3.12 or newer. Read the [genomic range contract](guides/ranges.md), [typed sequence Series contract](guides/sequences.md), [identifiers.org MCP guide](guides/identifiers.md), and [identifier-backed analysis guide](guides/artifacts.md) before using these APIs.
