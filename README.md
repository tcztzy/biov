# BioV

Biological tools, data and reproducible execution, available from Python, the
command line and MCP.

BioV's goal is a biology-focused tool manager: install or run a tool, reuse its
environment across tasks, and know which version ran and where its results live.
The model is `uv tool` and `uvx`, with biological contracts for reference versions,
coordinates, units and complete data.

## What works today

- Run scientific programs through declared, locked Pixi environments or explicit
  on-demand package sources
- Resolve biological identifiers to native provider files and reusable local
  artifacts, with Python and fsspec access
- Work with genomic intervals and explicitly typed DNA, RNA and protein sequences
- Run local managed Python analyses, retain complete outputs and records, and
  inspect bounded previews through Python, CLI or MCP
- Use the migrated `crisprprimer`, `crisprprimer-docker`, `biov-azimuth` and
  paired-read alignment interfaces

A unified installed-tool inventory, upgrades, uninstall and safe cache cleanup
are [planned](docs/guides/environments.md#planned-tool-lifecycle). They are not
implied by the current execution commands. In particular, `biov update` refreshes
the identifier registry, not installed tools.

BioV works without an agent or model. The caller chooses the scientific question,
method and interpretation; BioV handles repeatable data and execution rules.

## Install and run

The current package requires Python 3.12 or newer. From this checkout:

```sh
uv sync --locked
uv run biov --help
```

For Python use, import `biov` in that environment. To install the CLI from the
checkout into an isolated tool environment:

```sh
uv tool install .
biov --help
```

Set up the pinned Pixi manager and a declared scientific environment, then run
its native entry point:

```sh
biov setup
biov setup samtools
biov exec samtools --version
```

Declared environments install from their locks and can also be prepared on demand
by `biov exec`. Bundled scientific environments currently target Linux x86_64;
Pixi manager support on another platform does not establish support for every
scientific tool. BioV's wheel does not contain Pixi or the scientific packages.
See [scientific environments](docs/guides/environments.md) for platform support,
project manifests, offline manager archives and prerequisites.

`biov exec [SOURCE:]NAME ARGS...` supports `conda:`, `pypi:` and `npm:` sources.
A bare name selects a declared Pixi environment if present, otherwise a temporary
Pixi environment. PyPI and npm use `uv tool run` and `npx --yes`; those on-demand
sources do not use the project's lock. Missing managers produce an error rather
than silently changing sources. `biov pixi` passes native Pixi arguments through.

Run an ordinary script with `biov run analysis.py`; use
`biov run --executor lsf analysis.py` for an LSF submission receipt. A receipt is
not a completion result. Managed analysis is a separate, currently local path.

## Python data APIs

This example runs locally without downloading data:

```python
import biov
import pandas as pd

regions = biov.BioDataFrame({"seqid": ["chr1"], "start": [0], "end": [10]})
mask = biov.BioDataFrame({"seqid": ["chr1"], "start": [5], "end": [15]})
clipped = regions.intersect(mask)  # chr1: [5, 10)

dna = pd.Series(["ACGT", None], dtype="biov.dna")
reverse = dna.seq.reverse_complement()
gc = dna.seq.gc_fraction()
```

Interval APIs use 0-based, end-exclusive coordinates. They preserve row order and
duplicates under their documented rules. Downloaded files retain their native
coordinates; a download alone is not a coordinate conversion. Sequence types are
explicit, so ordinary strings are not guessed to be DNA or protein. Read the
[range](docs/guides/ranges.md) and [sequence](docs/guides/sequences.md) contracts.

Identifier-backed files are resolved on the host that uses them:

```python
import biov

protein_file = biov.path("uniprot://P42212")
with biov.open("uniprot://P42212", artifact="entry_json", mode="rt") as handle:
    metadata = handle.read()
```

These calls can download from the official provider on a cache miss. Use
`biov.artifact_capabilities()` to inspect supported namespaces and formats
without network access. RefSeq downloads require the official NCBI Datasets CLI
on the execution host. See [provider files](docs/guides/artifacts.md) for fsspec,
cache behavior and exact version handling.

## Complete results, bounded previews

`biov analyze REQUEST.json` runs a declared local analysis;
`biov inspect-analysis RECORD` reads its saved facts without resubmitting it.
The equivalent Python and MCP APIs share the same implementation.

Outputs and run records persist under `BIOV_ANALYSIS_ROOT`. Default responses
show at most two preview rows/records within a 32 KiB response cap; downstream
analysis consumes the complete saved file. Missing or changed outputs fail
explicitly. Execution success, scientific checks and biological interpretation
are separate facts.

MCP can read small registered outputs up to its 1 MiB resource limit. Larger
client downloads require configured existing HTTP(S) storage or explicit transfer;
an executor-local path is not automatically available to another client. See
[managed analysis](docs/guides/analysis.md) for the two-step example, security
boundary, tested clients/platforms and current remote-execution limits.

## MCP and agent plugins

Start the stdio server with `biov mcp`. It exposes identifier discovery, provider
records and managed analysis tools/resources. It does not require a wrapper tool
for each Python function. See the [MCP guide](docs/guides/identifiers.md).

The repository also supplies Codex and Claude Code plugins with shared skills.
[Plugin installation](docs/guides/plugins.md) is separate from installing BioV and
its scientific environments. Skills guide method selection and interpretation;
reusable computation stays in BioV or established scientific tools.

## Rust direction

The accepted target is a shared Rust core with a thin Python interface, plus CLI
and MCP access. The present implementation is still predominantly Python with
pandas, Biopython and RuRanges' Rust-backed kernels.

Polars is the preferred dataframe candidate. Rust-Bio, noodles and direct Rust
interval kernels will be evaluated against BioV's scientific contracts. Python
remains a first-class interface, but this actively developed package may change
names, signatures and return types, including replacing pandas/Biopython objects.
No legacy compatibility layer is required merely to preserve an old API.

This language choice is motivated by AI-assisted engineering and stronger
compile-time checks, independently of uv's language choice. Scientific correctness
still requires independent checks. See the [Rust migration plan](docs/guides/rust-migration.md)
and [SPEC](SPEC.md) for stages, breaking-change policy and native-distribution gates.

## Configuration and development

`BIOV_HOME` selects the data cache root; `BIOV_ANALYSIS_ROOT` selects persistent
analysis storage. Software environments and saved results have separate retention
rules. See [configuration](docs/guides/configuration.md) for TOML, environment
variables, fsspec settings and execution hosts.

```sh
uv run --locked pytest tests/ -q
uv run --locked mkdocs build --strict
uv run --locked prek run --all-files
uv build
```

The current source distribution includes package sources, resources, tests and
documentation. Plugin files are distributed through Git. Native Rust packaging
is planned, not part of this release. Existing score formulas and models are not
biologically validated merely by migration; see the
[CRISPR computation guide](docs/guides/crispr-computation.md).

[Documentation](https://tcztzy.github.io/biov/) · [MIT license](LICENSE)
