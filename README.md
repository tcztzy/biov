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
- Use Rust-native locked setup and literal execution for initial Samtools/GOATOOLS
  entry points, with [recorded local setup inspection](docs/guides/environments.md#native-rust-locked-tool-migration)
- Resolve biological identifiers to native provider files and reusable local
  artifacts, with Python and fsspec access
- Retain explicitly selected Hugging Face [model resource files](docs/guides/model-resources.md)
  at a full immutable Git commit through official `hf`, with portable
  checksums/README and offline complete-byte inspection; no model execution
- Work with genomic intervals and explicitly typed DNA, RNA and protein sequences
- Run local managed Python analyses, retain complete outputs and records, and
  inspect bounded previews through Python, CLI or MCP
- Open local CSV tables, query complete data with Rust Polars, and export/retrieve
  Arrow IPC through the Python-free [native MCP slice](docs/guides/rust-datasets.md);
  reopen saved native artifacts with their records after a server restart
- Generate typed GC-window datasets from one exact [prepared RefSeq
  FASTA](docs/guides/prepared-fasta.md#native-sequence-metric-tables), with explicit
  canonical and weighted-IUPAC policies; use the same native query/export/reopen workflow
- Use the migrated `crisprprimer`, `crisprprimer-docker`, `biov-azimuth` and
  paired-read alignment interfaces

`biov install`, `biov list` and `biov uninstall` are thin native adapters to
Pixi global for Samtools and uv tool for GOATOOLS on Linux x86_64. Their upstream
commands and environments live in dedicated directories under BioV's environment
root, separate from bundled locked scientific workflows. Listing shows the
managers' human-readable inventory in those roots. Upgrades and general cache
cleanup remain planned. `biov update` refreshes the identifier registry, not
installed tools.

BioV works without an agent or model. The caller chooses the scientific question,
method and interpretation; BioV handles repeatable data and execution rules.

## Install and run

The current package requires Python 3.12 or newer. Source builds now also require
Rust 1.89.0 (pinned in `rust-toolchain.toml`) and a C linker. Install Rust through
[the official rustup installer](https://rustup.rs/), then from this checkout:

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

The installed `biov` executable is Rust-owned. The wheel includes the existing
Python modules and native extension; Python-backed commands use the interpreter
paired with that installation. There is no separate `biov-rs` command.

With existing Pixi 0.81.0 and uv, install a supported scientific tool and expose
its upstream commands:

```sh
biov install samtools
biov install goatools
biov list
samtools --version
goatools find_enrichment --help
biov uninstall goatools
```

BioV prints PATH instructions for `environment_root/pixi-global/bin` and
`environment_root/uv-tools/bin` when needed. Apply them in your shell; BioV does
not edit startup files or user-global manager settings. Samtools is pinned to
1.24; GOATOOLS is pinned to 1.6.5 with statsmodels 0.14.6, using the Python
interpreter paired with the installed BioV package. The managers resolve remaining
dependencies; these installations do not consume the bundled scientific lock.

Upstream commands and environments remain usable after `uv tool uninstall biov`;
the GOATOOLS environment still requires its base Python installation to remain.
`biov uninstall NAME` asks the backend to remove that tool's environment and
exposed commands. Scientific data, saved outputs, caches and separate locked
workflow environments are unaffected. See the [scientific environments guide](docs/guides/environments.md)
for root selection, manager prerequisites and limits.

Bundled scientific environments currently target Linux x86_64. BioV's wheel
contains neither Pixi nor scientific packages. Existing manager bootstrap and
broader environment provisioning remain available explicitly through
`biov python setup`, including offline archives and project manifests. This
compatibility route does not populate the separate Pixi global or uv tool
inventory. See the scientific environments guide for prerequisites and scope.

`biov exec [SOURCE:]NAME ARGS...` supports `conda:`, `pypi:` and `npm:` sources.
A bare name selects a declared Pixi environment if present, otherwise a temporary
Pixi environment. PyPI and npm use `uv tool run` and `npx --yes`; those on-demand
sources do not use the project's lock. Missing managers produce an error rather
than silently changing sources. `biov pixi` passes native Pixi arguments through.

The native `biov tools exec` route provisions/reuses the bundled locked Samtools
or GOATOOLS environment without exposing an upstream user command or populating
the separate global-tool inventory. `--no-install` requires an existing native
provisioning receipt.
Use `biov tools inspect NAME` to inspect that provisioning state. The broader
`biov exec` route retains its existing Python-backed behavior.

`BIOV_ENVIRONMENT_ROOT` and `BIOV_PIXI_BIN` select native storage and the matching
Pixi; `BIOV_UV_BIN` or `--uv FILE` selects uv for the install/list/uninstall
adapters. Native routes do not read the Python TOML configuration. A standalone
`cargo install --locked --path crates/biov-cli` provides native routes, including
Samtools installation; GOATOOLS installation requires the complete wheel's paired
Python interpreter. Install the complete wheel for Python-backed capabilities.

Run an ordinary script with `biov run analysis.py`; use
`biov run --executor lsf analysis.py` for an LSF submission receipt. A receipt is
not a completion result. Managed analysis is a separate, currently local path.

## Selected model files

The Rust-native model route delegates acquisition to official `hf`. It prefers a
compatible installed client; otherwise an existing uv supplies
`huggingface-hub==2.1.1` with upstream SOCKS support in its on-demand environment.
Initial downloads support Linux x86_64. Select exact files and a full immutable
commit, using a new destination:

```sh
biov model download \
  --revision f171d7baecaf37b5da5a3616d8833b9969753535 \
  --local-dir ./tiny-bert-config \
  hf-internal-testing/tiny-random-bert config.json tokenizer_config.json
biov model inspect ./tiny-bert-config
```

These example commands select only configuration, not weights or a complete
model. Inspect and exact verified reuse are offline; the saved native files,
relative checksum record and standard-reader README remain usable after moving
the directory and removing BioV. No downloaded model code is executed, and
checksums do not establish upstream authenticity or scientific quality. See the
[model resource guide](docs/guides/model-resources.md) for client selection,
`--no-install`, failure retention, format bounds and current acceptance status.

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
and MCP access. Normalization, IUPAC reverse complement, validated sequence
lengths and weighted GC fractions execute in the shared Rust core through PyO3. The rest remains predominantly Python with pandas,
Biopython and RuRanges' Rust-backed kernels. The new batch API returns nullable
Python lists of strings or numeric results; see the [sequence contract](docs/guides/sequence-contract.md).

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
cargo test --workspace --locked
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
uv run --locked pytest tests/ -q
uv run --locked mkdocs build --strict
uv run --locked prek run --all-files
uv build
```

An advisory [AgentRC](https://github.com/microsoft/agentrc) report runs separately
from those engineering gates. Its prerelease CLI and transitive dependencies are
pinned in `.github/agentrc/package.json` and `package-lock.json`, outside BioV's
runtime dependencies. With Node.js 22+ and npm:

```sh
npm --prefix .github/agentrc ci --ignore-scripts --omit=optional --no-audit --no-fund
npm --prefix .github/agentrc test
npm --prefix .github/agentrc run --silent report > /tmp/biov-readiness.json
```

The local policy checks Rust/Python command and configuration presence, including
both lockfiles. It does not execute or certify builds, scientific tests or type
checking. The report may still suggest irrelevant Node/agent tooling. No maturity
level or pass-rate threshold gates CI; generated JSON is an uncommitted CI artifact.
Official `agentrc init` was inspected, not executed: its default instruction
generation uses Copilot and its other selections create editor/MCP settings.
The official-schema `agentrc.config.json` adds the Python source area alongside
automatically detected Cargo crates. Readiness needs no model login or calls.

The current source distribution includes package sources, resources, tests and
documentation, Cargo sources and the lockfile. Plugin files are distributed through
Git. Setuptools-rust builds the native extension and the Rust `biov` executable. Installing a BioV wheel does not compile
BioV, but dependencies may still need compilers where their wheels are unavailable.
This slice validates local Linux x86_64 builds, not a portable release matrix.
See [source builds and release gates](docs/guides/rust-migration.md#source-builds-and-native-validation).
Existing score formulas and models are not
biologically validated merely by migration; see the
[CRISPR computation guide](docs/guides/crispr-computation.md).

[Documentation](https://tcztzy.github.io/biov/) · [MIT license](LICENSE)
