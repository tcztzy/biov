BioV
====
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Ftcztzy%2Fbiov%2Fmain%2Fpyproject.toml)
![PyPI - Downloads](https://img.shields.io/pypi/dd/biov)

Next-generation development experience for computational molecular biology.

## Highlights

- **LLM first**: Designed for seamless integration with large language models, built for LLM workflows, and optimized for LLM-assisted development
- **Pydantic-powered**: Built-in validation and serialization for robust data handling
- **Pandas ecosystem**: Developer-friendly DataFrame operations with extended bioinformatics capabilities
- **RuRanges interval kernel**: Stable `BioDataFrame` range semantics over NumPy/Rust kernels
- **Typed sequences**: Explicit nullable DNA, RNA, and protein Series with a `.seq` API
- **Modern tooling**: Full type hints support and configuration through environment variables

## Coordination system
> [!IMPORTANT]
> BioV consistently uses BED-like, 0-based, end-exclusive `[start, end)` coordinates, regardless of input format (including GFF3 and VCF).

This design decision was made to (by Gemini 2.5 Pro Exp):
1. Direct Compatibility: It aligns seamlessly with Python slicing and the indexing conventions of most relevant programming languages.
2. Reduced Errors: Minimizes the risk of off-by-one errors, which are notoriously common when converting between 1-based/inclusive and 0-based/semi-open systems.
3. Simplicity: Length calculation (end - start) and handling adjacent/empty intervals are mathematically cleaner and more intuitive within a programming context.
4. Developer Familiarity: Most developers working with sequences in code are already accustomed to this paradigm.

## Requirements and interval engine

BioV requires Python 3.12 or newer. Genomic interval methods on `BioDataFrame` use BioV-owned pandas/NumPy adaptation around RuRanges' Rust kernels. PyRanges objects and conversion helpers are not part of the API.

Range operations group by chromosome and, when present on both operands, exact `+`/`-` strand. They preserve input order and duplicate rows. See the [genomic range contract](docs/guides/ranges.md) for overlap, intersection, subtraction, nearest-direction, empty-input, and coordinate details.

## Typed sequence Series

Importing BioV registers three explicit pandas extension dtypes. Ordinary string Series are never guessed to be biological sequences.

```python
import pandas as pd
import biov

dna = pd.Series(["ACGT", None], dtype="biov.dna")
dna.seq.reverse_complement()
dna.seq.gc_fraction()

protein = pd.Series(["ACDE"], dtype="biov.protein")
protein.seq.molecular_weight()
```

DNA/RNA reverse complement, weighted GC, and translation plus protein molecular weight, isoelectric point, and amino-acid composition use Biopython. See the [typed sequence contract](docs/guides/sequences.md) for alphabets, missing values, and error behavior.

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
