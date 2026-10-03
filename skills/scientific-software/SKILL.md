---
name: scientific-software
description: Find scientific software and run Python, R, or native tools through BioV. Use when choosing a biological analysis package, discovering available software, checking execution requirements, or running analysis locally or remotely.
---

# Scientific software and execution

Use the relevant task skill and the selected package's official documentation
when choosing software. Follow documentation for the installed version's API
and arguments; task guidance does not establish that a package is installed.

Before writing an unfamiliar call, read [native invocation reference](references/invocation.md)
to obtain import paths, signatures, defaults and command options from the selected
execution environment. Historical Biomni function names are not installed APIs.

Read the [execution guide](../../docs/guides/environments.md) for assigning
software to Pixi environments by manifest environment name, and declaring packages
for a new environment. Distinguish a declared environment, a successful
startup, and a completed analysis. Model weights and reference data are separate
inputs.

## Execute

The host must have the `biov` command installed. For a user-requested environment
setup, read the execution guide, run `biov setup` on the execution host, declare
the required packages in `[tool.pixi.*]` in `pyproject.toml`, generate its
`pixi.lock` with Pixi, and run `biov setup ENVIRONMENT`. Validate the native
command or import before analysis; environment creation alone is not scientific
validation. The shipped manifest uses conda-forge and Bioconda; do not add `defaults`
or Anaconda repository URLs. Preserve the manager and package license notices.
For local Python analysis with retained results, use MCP `run_analysis` and
`inspect_analysis` as described in the [managed analysis guide](../../docs/guides/analysis.md).
That path uses a declared locked Pixi environment and returns bounded previews
plus complete result references. Use the host agent's terminal tool for native
commands, R and the existing SSH/LSF interfaces; managed analysis does not
dispatch through SSH or LSF.

```sh
biov exec python analysis.py
biov exec r analysis.R
biov exec scvi analysis.py
biov exec samtools faidx reference.fa
biov exec pypi:ruff check analysis.py
```

Use native package APIs inside the script. Pass input and output paths explicitly,
retain native file formats, and inspect the actual exit status and output files.
Do not substitute a different program or fabricate a result after a failed run.

Deployment selects the execution host. A bare name defaults to conda and selects
a declared same-name environment, otherwise temporary Pixi execution. Declared
environments install from their lock and run their declared preparation task once
on demand.

Execution selects a same-name native Pixi task before a same-name
executable. The [shipped manifest](https://github.com/tcztzy/biov/blob/main/src/biov/assets/environments/pyproject.toml)
defines those entry tasks. Library tasks such as `scvi` and `deeppurpose` run
Python scripts; their `--help` is Python help, not a library CLI or an import
check. `r` invokes Rscript. Select another program in a suite with native
`biov pixi run` arguments.

`conda:NAME` uses that declared environment if present, otherwise temporary Pixi execution;
`pypi:NAME` and `npm:NAME` use uv and npx. Temporary sources resolve without the
project lock. `--no-install` skips preparation and installation only for declared
environments and rejects temporary sources. Run commands already on the host PATH
directly in your shell. Missing managers, executables, weights or reference data
fail without substituting another source or program.

Read the [execution guide](../../docs/guides/environments.md) for working
directories, shared input/output paths and environment behavior. Pixi and SSH
configuration belong to the deployment operator; the analysis invocation remains
`biov exec [SOURCE:]NAME ARGS...`, with BioV options before the coordinate.

For input discovery, read [biological-data](../biological-data/SKILL.md).
Use `biov.path`, `biov.open`, or fsspec inside the execution environment to obtain
files. For scientific decisions and QC, load the relevant analysis skill from
the host's skill list. Read only the instructions and reference documents
needed for the task; no separate LLM resource-selection step is required.
