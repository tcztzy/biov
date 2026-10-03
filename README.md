BioV
====
![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Ftcztzy%2Fbiov%2Fmain%2Fpyproject.toml)
![PyPI - Downloads](https://img.shields.io/pypi/dd/biov)

BioV builds high-performance bioinformatics infrastructure for AI agents. Its
Python APIs provide genomic interval and sequence operations and access to
identifier-backed data; its MCP server exposes identifier discovery, provider
records and [local managed Python analysis](docs/guides/analysis.md) with saved
results; its CLI also runs scripts directly or submits them to LSF.

BioV provides reusable capabilities that work without an agent. Task-specific
tool orchestration, decisions about the next analysis step, and biological
interpretation belong to the agent or workflow project using BioV.

## Highlights

- **Agent-facing interfaces**: MCP tools for identifier discovery, provider records and local managed analysis, with bounded previews and access to complete saved results
- **Pydantic-powered**: Built-in validation and serialization for robust data handling
- **Pandas ecosystem**: Developer-friendly DataFrame operations with extended bioinformatics capabilities
- **RuRanges interval kernel**: Stable `BioDataFrame` range semantics over NumPy/Rust kernels
- **Typed sequences**: Explicit nullable DNA, RNA, and protein Series with a `.seq` API
- **Persistent identifiers**: Resolve identifiers.org Compact Identifiers through MCP resources and prompt parsing
- **Portable analysis**: Resolve persistent IDs as executor-local `PathLike` artifacts, then use ordinary Biopython or command-line tools
- **Modern tooling**: Full type hints support and configuration through environment variables

## Deterministic CRISPR computation

BioV distributes the existing `crisprprimer` Python package and CLI previously
owned by GEEPilot, together with `crisprprimer-docker` and `biov-azimuth`. These
provide repeatable computation and native report parsing. GEEPilot task skills
retain question selection, safety routing, method choice and interpretation.
The migration preserves existing score formulas and presets; it does not establish
their biological validity. See the [computation guide](docs/guides/crispr-computation.md)
for interfaces, environments and verification limits.

## Coordinate system
> [!IMPORTANT]
> BioDataFrame interval APIs use BED-like, 0-based, end-exclusive `[start, end)` coordinates. Downloaded provider files retain their native coordinate conventions.

The interval convention matches Python slicing, with length `end - start`.
Use the format's parser to convert coordinates before applying interval APIs;
reading or downloading a raw file does not itself normalize its coordinates.

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

`biov mcp` exposes provider data through
`refseq.gcf://GCF_000001030.2`, `uniprot://P42212`,
`pubmed://22140103`, `clinvar://65533`, `dbsnp://rs121909098`, and
`geo://GSE1000`. Registry metadata and
resolver responses use `identifiers://<registry>` and
`identifiers://<registry>:<id>`. The `parse_identifiers` tool recognizes
Compact Identifiers, identifiers.org URLs, these resource URIs, and explicitly
allowlisted unambiguous bare IDs such as `GCF_000001030.2`. The
`resolve_identifiers` tool embeds the same resource content for MCP clients
that cannot call `resources/read`. Refresh the raw registry response with
`biov update`. See the
[identifiers.org MCP guide](docs/guides/identifiers.md) for exact resource
contents, host configuration, and error behavior.

The file namespaces and available representations are described by
`biov.artifact_capabilities()` and the [file provider guide](docs/guides/artifacts.md).
MCP resources for these providers return JSON describing how to open their files;
RefSeq and UniProt retain their native metadata resources. Actual data is read
through `biov.path`, `biov.open`, or fsspec inside the analysis environment.

Database searches and biological analysis belong to the task
[skills](skills/) and their upstream libraries or services. BioV no longer
publishes API catalogs or a generic `query_database` tool. The
[biological-data skill](skills/biological-data/SKILL.md) directs explicit queries
to official API documentation. See the [file provider guide](docs/guides/artifacts.md)
for supported file formats and limitations.

## Identifier-backed analysis

Keep persistent IDs in generated code and resolve them to paths inside the selected
execution environment. BioV owns the data boundary; established libraries own
the computation:

```python
import biov
from Bio import SeqIO

records = SeqIO.parse(biov.path("GCF_000006945.2"), "fasta")
```

BioV registers every namespace in its artifact capability manifest with fsspec:

```python
import fsspec

with fsspec.open("uniprot://P42212", "rt") as sequence_file:
    fasta_text = sequence_file.read()
```

These reads reuse `biov.path` downloads and cache files. Existing readers work
with the same URIs: `biov.read_fasta("uniprot://P42212")` still returns sequence
records, while `biov.read_gff3(uri, storage_options={"artifact": "annotation_gff3"})`
returns a `BioDataFrame` for a RefSeq URI.

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

## Agent plugins

This repository also provides a `biov` plugin marketplace for Codex and Claude
Code. Both hosts share the same task skills, reference files and MCP server.
The `scientific-software` skill explains how to use configured Pixi environments
and invoke scientific software; skills and references are loaded on demand.
See the [plugin installation guide](docs/guides/plugins.md) for runtime
prerequisites, local/GitHub installation and validation. Plugin installation is
separate from installing the Python package and scientific software.

## Environments

Development uses `pyproject.toml`, `uv.lock` and
[uv](https://docs.astral.sh/uv/getting-started/installation/):

```sh
uv sync --locked
uv run --locked pytest tests/ -q
uv run --locked prek run --all-files
uv build
```

The source distribution includes the Python package and its resources, Python
tests, documentation and documentation build files. Skills and plugin files are
distributed through Git; their tests run from the repository.

The uv development environment contains BioV's declared dependencies and development
tools. Scientific scripts and native programs run through
`biov exec [SOURCE:]NAME ARGS...`, for example `biov exec python analysis.py`.
A bare name defaults to conda: it first selects a declared Pixi environment,
otherwise a temporary Pixi environment. Declared environments install from their
lock and run their declared `prepare` task once before execution. `conda:NAME` uses the same
locked environment if declared, otherwise a temporary Pixi environment;
`pypi:NAME` uses `uv tool run`, and `npm:NAME` uses `npx --yes`. These temporary
sources resolve packages without the project's lock. Missing uv or npx produces
an error with installation instructions; BioV does not switch sources.

In a declared environment, a same-name native Pixi task takes precedence over
the same-name executable. Entry tasks live in the shipped manifest:
`biov exec r analysis.R` invokes Rscript, and `biov exec scvi analysis.py` runs a Python script
in the scvi environment. `biov exec vina-meeko ARGS...` invokes Vina; use
`biov pixi run` to select other programs in that suite.

Arguments pass through without requiring `--`. Put BioV's `--cwd` and
`--no-install` before the coordinate; `--no-install` skips installation and
preparation for declared environments and rejects temporary sources. BioV reads application and SSH
configuration on the execution host. Native Pixi commands remain available
through `biov pixi`. Run commands already on the host PATH directly in your shell.
Run `biov setup` to make a Pixi manager available on a new execution
host, then `biov setup ENVIRONMENT` to install a locked scientific environment.
A `BIOV_PIXI_BIN` path, then a matching `pixi` on `PATH`, is reused; only the
managed copy under `BIOV_ENVIRONMENT_ROOT` is downloaded. Pixi is an optional
runtime prerequisite for Pixi-backed commands and does not manage BioV
development. Scientific requirements live in `[tool.pixi.*]` in
`src/biov/assets/environments/pyproject.toml`, with resolved
versions in its adjacent `pixi.lock`. Scientific environments target Linux x86_64.
`biov setup --all` installs all declared environments without an agent. See the
[software guide](docs/guides/environments.md).

## Configuration

Set environment variables before starting BioV:

- `BIOV_HOME`: Path to custom cache directory (default: platform-specific cache dir)
- `BIOV_CACHE_HTTP`: Enable/disable HTTP caching (default: True)

The cache directory is determined by:
1. `BIOV_HOME` if set
2. `XDG_CACHE_HOME/biov` if XDG_CACHE_HOME is set
3. Platform-specific cache directory otherwise

Existing fsspec `filecache` configuration and explicit `cache_storage` options
take precedence over BioV's default for ordinary cached URLs. Local data paths
are configuration, not packaged metadata. See the
[configuration guide](docs/guides/configuration.md) for cache examples.

## Executables

- biov

Native scientific programs run through `biov exec` using package coordinates;
BioV does not install a wrapper per program.

## Supported formats

- [x] GFF3
- [x] PSL
- [x] FASTA
- [ ] BED
- [ ] VCF
