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
- **Persistent identifiers**: Resolve identifiers.org Compact Identifiers through MCP resources and prompt parsing
- **Portable analysis**: Resolve persistent IDs as executor-local `PathLike` artifacts, then use ordinary Biopython or command-line tools
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

## Identifiers.org MCP server

`biov mcp` exposes canonical data through
`refseq.gcf://GCF_000001030.2` and `uniprot://P42212`. Registry metadata and
resolver responses use `identifiers://<registry>` and
`identifiers://<registry>:<id>`. The `parse_identifiers` tool recognizes
Compact Identifiers, identifiers.org URLs, these resource URIs, and explicitly
allowlisted unambiguous bare IDs such as `GCF_000001030.2`. The
`resolve_identifiers` tool embeds the same resource content for MCP clients
that cannot call `resources/read`. Refresh the raw registry response with
`biov update-identifiers-registry`. See the
[identifiers.org MCP guide](docs/guides/identifiers.md) for exact resource
contents, host configuration, and error behavior.

## Identifier-backed analysis

Keep persistent IDs in generated code and resolve them to paths inside the selected
execution environment. BioV owns the data boundary; established libraries own
the computation:

```python
import biov
from Bio import SeqIO

records = SeqIO.parse(biov.path("GCF_000006945.2"), "fasta")
```

The `refseq.gcf × genome_fasta` provider invokes the official NCBI
`datasets download genome accession` command and preserves the complete
extracted package. UniProt JSON and FASTA representations are fetched from
their official `/uniprotkb/<accession>.json|.fasta` endpoints and cached
independently without rewriting them. The artifact kind defaults from the namespace, so
`biov.path("uniprot://P42212")` returns `P42212.fasta`; request `entry_json`
for the complete metadata and database cross-references. Cache fills are atomic
and valid cache hits skip the corresponding downloader. Install the
[NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
on execution hosts that resolve RefSeq genome paths.

Run the complete ordinary Python script locally with `biov run analysis.py`, or
submit the same script with `biov run --executor lsf analysis.py`. LSF returns a
job-ID submission receipt, not a false completion result. See the
[identifier-backed analysis guide](docs/guides/artifacts.md) for Biopython code,
cache behavior, and HPC deployment requirements.

## Environments

BioV can be configured through environment variables (prefixed with `BIOV_`) or a `.env` file:

- `BIOV_HOME`: Path to custom cache directory (default: platform-specific cache dir)
- `BIOV_CACHE_HTTP`: Enable/disable HTTP caching (default: True)

The cache directory is determined by:
1. `BIOV_HOME` if set
2. `XDG_CACHE_HOME/biov` if XDG_CACHE_HOME is set
3. Platform-specific cache directory otherwise

## Executables

- biov
- blat

## Supported formats

- [x] GFF3
- [x] PSL
- [x] FASTA
- [ ] BED
- [ ] VCF
