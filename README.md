BioV
====
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Ftcztzy%2Fbiov%2Fmain%2Fpyproject.toml)
![PyPI - Downloads](https://img.shields.io/pypi/dd/biov)

Next-generation development experience for computational molecular biology.

## Highlights

- **LLM first**: Designed for seamless integration with large language models, built for LLM workflows, and optimized for LLM-assisted development
- **Pydantic-powered**: Built-in validation and serialization for robust data handling
- **Pandas ecosystem**: Developer-friendly DataFrame operations with extended bioinformatics capabilities
- **Modern tooling**: Full type hints support and configuration through environment variables

## Coordination system
> [!IMPORTANT]
> BioV consistently uses BED-like coordinates with 0-based start positions and 1-based end positions, regardless of input format (including GFF3 and VCF).

This design decision was made to (by Gemini 2.5 Pro Exp):
1. Direct Compatibility: It aligns seamlessly with Python slicing and the indexing conventions of most relevant programming languages.
2. Reduced Errors: Minimizes the risk of off-by-one errors, which are notoriously common when converting between 1-based/inclusive and 0-based/semi-open systems.
3. Simplicity: Length calculation (end - start) and handling adjacent/empty intervals are mathematically cleaner and more intuitive within a programming context.
4. Developer Familiarity: Most developers working with sequences in code are already accustomed to this paradigm.

## Environments

BioV can be configured through environment variables (prefixed with `BIOV_`) or a `.env` file:

- `BIOV_HOME`: Path to custom cache directory (default: platform-specific cache dir)
- `BIOV_CACHE_HTTP`: Enable/disable HTTP caching (default: True)

The cache directory is determined by:
1. `BIOV_HOME` if set
2. `XDG_CACHE_HOME/biov` if XDG_CACHE_HOME is set
3. Platform-specific cache directory otherwise

## Executables

- blat

## Supported formats

- [x] GFF3
- [x] PSL
- [x] FASTA
- [ ] BED
- [ ] VCF
