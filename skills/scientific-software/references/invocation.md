# Native invocation reference

Use the task skill to choose the scientific method, then inspect only the native
function or command needed for that task. Run inspection through the same
package coordinate or native Pixi environment as the analysis: the host Python
and a specialist environment can contain different packages and versions. Retain the returned
help in context while constructing the call.

The examples below assume a project declaring the named environments and entry
tasks from the [shipped manifest](https://github.com/tcztzy/biov/blob/main/src/biov/assets/environments/pyproject.toml).
Select the same manifest for setup and native Pixi execution:

```sh
export BIOV_ENVIRONMENT_MANIFEST=/path/to/project/pyproject.toml
```

## Python: imports, parameters and defaults

Use BioV's application interpreter for its own APIs. Save this as
`inspect_biov.py` and run `biov run inspect_biov.py` on the execution host:

```python
import json
import biov
from importlib.metadata import version

print(version("biov"))
help(biov.path)
help(biov.open)
print(json.dumps(biov.artifact_capabilities(), indent=2))
```

Scientific environments contain the declared analysis packages, not another
BioV application installation. Prepare files with `biov run prepare.py`, then
pass their paths to `biov exec TOOL ARGS...` on a host that can access those files.

For a selected analysis library, inspect its actual callable. For example, in
the declared single-cell environment of an explicit project:

```sh
biov exec scvi -c 'import scanpy; from importlib.metadata import version; print(version("scanpy")); help(scanpy.pp.calculate_qc_metrics)'
```

`help` supplies the signature and docstring when available, including defaults,
input requirements and return values. Use the official documentation matching
the installed version for compiled functions or APIs with incomplete local help.
The distribution name used by `version` can differ from the import name
(for example, `biopython` and `Bio`). Inspect the selected model's documented
checkpoint and preprocessing requirements separately; a signature cannot supply them.

Save the analysis as a script and run `biov exec scvi analysis.py` in the selected
library environment, or `biov exec python analysis.py` in the general scientific
Python environment. Library entry tasks invoke Python, so `--help` alone does not
verify a library import or its API. Scientific tools receive files prepared in
the BioV application environment;
do not assume that BioV is importable in a specialist runtime.

## R: package functions and help

Use Rscript in the declared `r` environment. Inspect the selected package without
starting an interactive help browser:

```sh
biov exec r -e 'options(pager="cat"); print(packageVersion("stats")); print(args(stats::lm)); print(help("lm", package="stats", help_type="text"))'
```

For DESeq2, substitute `DESeq2` for the package and `DESeq` for the function
and help topic. Then run the saved script with
`biov exec r differential_expression.R`.
For S4 generics, consult the documented method for the actual input class.

## Native commands and suites

```sh
biov exec samtools --version
biov exec samtools faidx --help
biov setup hmmer
biov pixi run --manifest-path "$BIOV_ENVIRONMENT_MANIFEST" --environment hmmer -- hmmscan -h
```

Use the program's documented help/version option; these flags are not universal.
Some programs return a nonzero status when displaying usage, so inspect the
output as well as the status. A suite environment can contain several executable
names; use native Pixi arguments to select one. Do not assume an environment name
is a Python or R executable. A same-name Pixi task takes precedence; without one,
BioV runs the same-name real executable. Consult the manifest for the entry
command. For example, `biov exec vina-meeko --help` uses Vina; select Meeko with
`biov pixi run --manifest-path "$BIOV_ENVIRONMENT_MANIFEST" --environment vina-meeko -- mk_prepare_ligand.py --help`.

Run help and analysis in the same configured working directory. For a remote
execution host, scripts and inputs must already exist at paths visible there;
`biov exec` does not transfer local files. Use
`biov exec --cwd /shared/project python analysis.py`
when an explicit execution directory is required.

## Service APIs

For database searches, use the service links in
[biological-data](../../biological-data/references/services.md) to obtain the
current official schema or SDK documentation. Select the endpoint, required
fields and pagination from that source. Identifier-file retrieval instead uses
the installed BioV manifest above. Neither route requires a legacy Biomni wrapper.
