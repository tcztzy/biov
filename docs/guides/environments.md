# Scientific environments

## Install, list and remove native tools

The single public `biov` executable owns routing in Rust. On Linux x86_64,
`install`, `list` and `uninstall` are thin adapters to established upstream
managers: Pixi 0.81.0 global manages Samtools 1.24, while uv tool manages
GOATOOLS 1.6.5 with statsmodels 0.14.6. Installation/removal and nonempty inventory
need the corresponding existing manager; these commands do not bootstrap managers
or switch backends silently.

```bash
biov install samtools
biov install goatools
biov list
samtools --version
goatools find_enrichment --help
biov uninstall goatools
```

Options precede the tool name. All three commands accept `--environment-root DIR`,
`--pixi FILE` and `--uv FILE`; `BIOV_PIXI_BIN` and `BIOV_UV_BIN` select matching
manager executables. The Samtools route requires Pixi 0.81.0. The GOATOOLS route
requires uv and the Python interpreter paired with the installed BioV package,
so use the complete wheel rather than a standalone Cargo binary to install it.
It does not select an ambient PATH Python or download another interpreter.
Supported paired deployments are `uv tool install biov` and installing the complete
wheel with pip inside an ordinary virtual environment. Ad-hoc `--user`,
`--prefix` or `--target` layouts and relocated standalone binaries do not guarantee
a sibling Python; GOATOOLS installation and Python-backed routes fail explicitly
when it is absent. Pure native routes do not require that interpreter.

uv 0.12.9 passed packaging CI and uv 0.12.19 passed the local real-tool gates.
Other versions are not promised; an unsupported option produces the upstream
error and a BioV backend-failure diagnostic, without another backend fallback.

`--environment-root DIR` overrides `BIOV_ENVIRONMENT_ROOT`; otherwise the native
root is `$XDG_DATA_HOME/biov/environments` or
`$HOME/.local/share/biov/environments`. Within that root, the upstream-owned layout
is separate from ordinary user-global manager installations and the bundled
content-keyed scientific workspaces:

- Pixi global: `pixi-global/envs`, commands in `pixi-global/bin`, native manifest
  at `pixi-global/manifests/pixi-global.toml`, and cache in `pixi-global/cache`
- uv tool: environments/receipts in `uv-tools/tools`, commands in `uv-tools/bin`,
  cache in `uv-tools/cache`, and managed-Python directory in `uv-tools/python`

The Pixi adapter explicitly selects its owned native global manifest. It does
not fall back to the user's XDG/global manifest. uv's tool/bin/cache/Python
locations are similarly isolated. BioV does not edit shell profiles or user-global
manager configuration. There is no `--bin-dir` option or `BIOV_BIN_DIR` override;
use the reported dedicated bins. Installation prints shell-specific PATH
instructions when needed. Run them in the shell/process that will use the tools;
a child process cannot alter its parent shell's PATH.

Samtools installation delegates to `pixi global install`, selects the `samtools`
environment, exposes `samtools`, uses Linux-64 packages from conda-forge and
bioconda, and disables desktop shortcuts. Pixi owns its native trampoline and
activation. GOATOOLS installation delegates to `uv tool install goatools==1.6.5
--with statsmodels==0.14.6`, with the paired interpreter, no Python downloads and
no project configuration. uv exposes GOATOOLS' upstream console scripts, including
`goatools` and `find_enrichment.py`; BioV does not invent another wrapper.

Isolation here concerns installation directories and owned inventory, not a
sanitized process environment or isolated package index. uv configuration-file
discovery is disabled, but manager environment settings remain inherited,
including `UV_INDEX_URL`, `UV_DEFAULT_INDEX`, `UV_INDEX`, `UV_FIND_LINKS`,
`UV_OFFLINE` and applicable `PIXI_OFFLINE`/`PIXI_FROZEN` settings. They may alter
resolution sources or prevent a new install. BioV does not scrub proxy/authentication
or other backend settings; explicit tool/bin/cache roots and paired `--python`
remain selected. No full transitive-lock or public-index-only guarantee is implied.

These routes pin the primary tool versions and the declared statsmodels
compatibility dependency. Remaining dependencies are resolved by the upstream
manager. They do not consume the shipped scientific Pixi lock or promise a
hash-bound, fully locked transitive package set. Repeated install/reinstall uses
the backend's native behavior; a later resolution may change transitive packages.
Use the separate [bundled locked workflow](#native-rust-locked-tool-migration)
when that exact manifest/lock contract is required.

`biov list` invokes the managers' own human-readable list operations only for
existing dedicated roots. It retains their complete output under path-labeled
`Pixi global (...)` and `uv tool (...)` sections. Each manager's stdout is bounded
to 1 MiB; oversized, unreadable or non-UTF-8 output fails rather than becoming a
truncated inventory. It does not parse uv's text into a custom schema, emit JSON
or claim an independent readiness/integrity audit.
A first-time empty root reports that no native tools are installed without
requiring a manager or creating the root. The catalog,
execution-only workspaces, ordinary external uv/Pixi installations and biological
data caches are outside this inventory.

Upstream entrypoints and their environments are independent of the management
BioV wheel: `uv tool uninstall biov` does not remove Samtools or GOATOOLS installed
in these roots. There is no copied BioV runner, custom installed registry or
publication journal. The GOATOOLS virtual environment still needs its base Python
installation; deleting that interpreter can break it. These are software
installations, not relocatable biological data bundles.

Management install/list/uninstall requires the corresponding external manager to
remain available through `--pixi`, `--uv`, their BioV environment selectors or PATH.
Removing uv/Pixi itself prevents later management operations until it is restored;
it is not the same as removing the BioV management wheel. Existing upstream
entrypoints do not invoke the BioV executable.

`biov uninstall NAME` delegates to the selected backend, removing that owned
user-tool environment and its exposed commands. It does not delete the other
tool, package caches, scientific data, model weights, saved outputs or separate
bundled locked-workflow environments. General purge/cache cleanup and an explicit
upgrade interface are not implemented by these commands. BioV rejects redirected
backend locations and bounded foreign/modified-command collisions before mutation;
the upstream managers still own entrypoints and installation/removal semantics.
Native failure diagnostics propagate without force overwrites or a parallel
repair journal. These conservative guards require trusted roots; they do not
authenticate packages or defend against hostile coordinated metadata/file changes.

`biov tools exec NAME ...` remains a separate native locked
provisioning/execution path for the two supported tools; it does not expose user
commands or populate the global-tool inventory. `biov tools inspect NAME` inspects
its provisioning receipt rather than the upstream inventory. Other existing
capabilities still use the explicit Python bridge, including bare `biov mcp`,
`analyze`, `inspect-analysis`, `exec`, `run`, `pixi` and registry `update`.
Native datasets use `biov mcp-native`. Python-only setup remains available as
`biov python setup`; the rest of this guide labels that broader compatibility
behavior separately.

Upstream references: [Pixi global install](https://pixi.prefix.dev/latest/reference/cli/pixi/global/install/),
[Pixi global list](https://pixi.prefix.dev/latest/reference/cli/pixi/global/list/),
[Pixi global uninstall](https://pixi.prefix.dev/latest/reference/cli/pixi/global/uninstall/),
[uv tools](https://docs.astral.sh/uv/concepts/tools/) and
[uv directory controls](https://docs.astral.sh/uv/reference/environment/).

## Existing Python environment compatibility


BioV runs in the Python installation that already hosts the application. Pixi
manages separate scientific environments; BioV does not install a second copy
of itself. Environment definitions use native `[tool.pixi.*]` tables in
`src/biov/assets/environments/pyproject.toml`, with resolved versions in the
adjacent `pixi.lock`. The root project and development tools use uv.

## New execution host

```bash
biov python setup
```

This prefers a Pixi the host already has: the executable named by
`BIOV_PIXI_BIN`, then a `pixi` on `PATH` whose reported version matches the
pinned release, and only then the copy BioV manages in
`BIOV_ENVIRONMENT_ROOT/pixi-<version>`. Only the managed copy is downloaded, and
its platform-specific SHA-256 is verified before first execution. A `PATH` Pixi
with a different version is reported and skipped rather than used silently. The
default root is the platform user data directory's
`biov/environments`. No administrator privileges or shell profile changes are
needed. Supported manager platforms are macOS arm64/x86_64, Linux aarch64/x86_64,
and Windows arm64/x86_64; individual scientific packages have their own platform
requirements. Linux uses the upstream musl binary.

A distributor can supply the unmodified official release archive alongside its
installer and initialize without downloading the manager:

```bash
biov python setup --archive /path/to/pixi-aarch64-apple-darwin.tar.gz
```

The wheel includes Pixi's BSD-3-Clause notice and the environment manifest and
lock file. It does not include the Pixi binary or scientific packages. Offline
manager installation does not make uncached scientific packages available.

## Declare and lock scientific dependencies

The repository declares separate environments for the scientific tools used by
its skills. All bundled scientific environments currently target `linux-64`; the Pixi
manager itself supports the other platforms listed above. Each environment uses
`no-default-feature = true`, keeping BioV's application/development dependencies
out of specialist runtimes. Preinstall them explicitly, or let `biov exec`
install the selected declared environment on demand:

```bash
biov python setup samtools
biov python setup --all
```

`biov pixi` runs the resolved manager with arguments passed through unchanged, so
any Pixi command is available without BioV-specific configuration; add
`--manifest-path` to target the bundled workspace or a project of your own.
There is no LLM, agent plugin, or Docker requirement in this installation path.

The [shipped manifest](https://github.com/tcztzy/biov/blob/main/src/biov/assets/environments/pyproject.toml) is the
source for environment names, dependencies, entry tasks and preparation tasks;
its adjacent lock records resolved versions.

Sommer and its source-only R dependencies use pinned, checksum-verified CRAN
archives and native `R CMD INSTALL`. UCE and DiffDock use pinned upstream source
revisions. AutoSite uses a
checksum-verified upstream installer. Their declared `prepare` tasks are run by
`biov python setup` after Pixi installs dependencies. They are not replaced by similarly
named third-party packages. Model weights, databases and reference atlases remain
separate user-selected inputs. ChatNT's remote model code is not executed by
setup. License-restricted services such as NetOGlyc and IUPred are not represented
as freely redistributable packages.

DiffDock's inference environment follows the
[official environment at its pinned source revision](https://github.com/gcorso/DiffDock/blob/85c49b60d3e0b0182a59ee43a34a6d7036981284/environment.yml),
including ESMFold, OpenFold and the compiler/CUDA dependencies needed to build
OpenFold. The shipped manifest and lock contain the corresponding declarations
and resolved packages.

For example, add these tables to a project's `pyproject.toml` (and select the
platforms that its packages support):

```toml
[tool.pixi.workspace]
channels = ["https://conda.anaconda.org/conda-forge", "https://conda.anaconda.org/bioconda"]
platforms = ["linux-64"]

[tool.pixi.feature.science.dependencies]
python = "3.12.*"
numpy = "*"

[tool.pixi.environments]
python = { features = ["science"], no-default-feature = true }
```

For development, regenerate the lock through BioV itself and commit both files:

```bash
BIOV_ENVIRONMENT_MANIFEST=/path/to/pyproject.toml biov python setup --update-lock
```

This also installs Pixi if absent. Lock updates require an explicit project
manifest; ordinary setup uses the shipped lock without changing it. Use Pixi's
native features/environments for different software suites; BioV has no parallel
package specification format.

Deployment settings use Pixi's native controls. `PIXI_CACHE_DIR` selects a local
cache when the home directory is on a slow network filesystem. For slow source
builds during lock updates, `UV_LOCK_TIMEOUT` sets the cache-lock wait in seconds.
An operator can configure community-channel mirrors in the selected project's
`.pixi/config.toml`, for example:

```toml
[mirrors]
"https://conda.anaconda.org/conda-forge" = ["https://prefix.dev/conda-forge"]
"https://conda.anaconda.org/bioconda" = ["https://prefix.dev/bioconda"]
```

## Select the project and command

Machine settings remain in the platform user configuration directory's
`biov/config.toml`, or the file selected by `BIOV_CONFIG`:

```toml
environment_manifest = "/path/to/project/pyproject.toml"
```

```bash
biov python setup python
biov exec python -c 'import numpy; print(numpy.__version__)'
```

`setup python` installs the environment named in `[tool.pixi.environments]`.
A missing or stale lock fails; setup does not silently resolve new versions.
Repeated setup reconciles the environment with the lock. Update the manifest
and lock explicitly before requesting an environment update.

Omit `environment_manifest` to use the manifest and lock shipped in the wheel.
Setup copies them into a writable, content-addressed workspace beneath
`environment_root/workspaces`; it does not write into site-packages. Each
manifest/lock revision gets its own workspace. Editable/source installations
use the repository's same two files as the copy source. Explicit project
manifests are used in place, with environments under the project's `.pixi`
directory by default.

`biov exec NAME ARGS...` selects a declared environment named exactly like the
coordinate, installs that environment from the lock when it is missing
(`pixi install --locked`), runs the environment's declared `prepare` task once if
it has one, and then selects its same-name Pixi task. The completion record is
stored inside the environment it prepared, so cleaning or reinstalling that
prefix prepares it again instead of treating it as ready. If the entry task is
absent, it runs the same-name executable in the environment. Execution uses
`pixi run --frozen`; it never re-solves dependencies or rewrites the lock.
`biov exec --no-install NAME ARGS...` skips provisioning and preparation, and
runs with `pixi run --as-is`, using the installed dependencies without changing
the lock. A name that matches no declared environment uses temporary conda
execution, described below.
Preparation is resolved for the current platform: a feature task overrides the
workspace default, and `[tool.pixi.target.<platform>.tasks]` or
`[tool.pixi.feature.<name>.target.<platform>.tasks]` overrides the unscoped task
in the same table.
There are no BioV command mappings or default-runtime settings.

Entry tasks belong to the native Pixi feature, alongside dependencies and
`prepare`. The shipped `r` task invokes Rscript; library environments such as
`scvi` and `deeppurpose` invoke Python with the supplied script:

```bash
biov exec r analysis.R
biov exec scvi analysis.py
biov exec deeppurpose predictions.py
```

These Python tasks are script entry points, not package-specific CLIs. Consult
the shipped manifest for each entry command. A project can declare an entry
task using native Pixi syntax, for example:

```toml
[tool.pixi.feature.r.tasks]
r = 'cd "$INIT_CWD" && Rscript'
```

Pixi runs tasks from the project root. The shipped entry tasks use Pixi's
`INIT_CWD` to run in the caller's working directory, including `biov exec --cwd`.
Entries for prepared sources use `PIXI_PROJECT_ROOT` to locate those sources.
Preparation continues to run in the project workspace.

The DiffDock entry supplies the default inference configuration from its prepared
source directory. Use the native `--config /path/to/config.yaml` argument to
override it; relative input and output paths still use the caller's directory.

For other executables in a suite, use native Pixi arguments. The shipped
`vina-meeko` entry invokes Vina; select Meeko's ligand preparation program
explicitly, using the same manifest for setup and execution:

```bash
biov exec vina-meeko --help
biov python setup vina-meeko
biov pixi run --manifest-path /path/to/project/pyproject.toml --environment vina-meeko -- mk_prepare_ligand.py --help
```

## Package sources

`biov exec [SOURCE:]NAME ARGS...` accepts `conda:`, `pypi:` and `npm:`.
The default source for bare names is conda:

```bash
biov exec conda:jq --version
biov exec pypi:ruff --version
biov exec npm:prettier --version
```

Both bare `NAME` and `conda:NAME` use a declared same-name environment and its lock
when present.
Otherwise BioV delegates to `pixi exec -s NAME -- NAME ARGS...`, which runs
outside a workspace in a temporary environment without the project's manifest
or lock. `pypi:NAME` delegates to `uv tool run NAME ARGS...`; `npm:NAME` delegates
to `npx --yes NAME ARGS...`. These sources resolve packages through their native
managers and caches. Clean temporary Pixi environments with
`biov pixi clean cache --exec`.

Missing uv or npx fails with an error naming the missing manager and its
installation instructions. BioV never silently selects another source. Only
these known prefixes identify sources; other colon-containing names remain
ordinary package coordinates in the default conda source. Commands already on
the host PATH can be run directly in your shell.

Put BioV's `--cwd DIR` and `--no-install` before the coordinate. Everything after
it is a native argument, including options with those same names. An initial `--`
separator before native arguments is optional. `--no-install` works only for
declared environments; it rejects temporary sources before execution because
their managers may create environments or install packages.

Native arguments, working directory, stdio, and exit status are preserved. Caller
cache variables remain inherited; BioV supplies `BIOV_HOME` and
`BIOV_CACHE_HTTP` for every source. Configure package environments with their
native manager controls. Missing prerequisites fail without switching sources.

Use `biov run prepare.py` on the execution host for data preparation using
BioV APIs, then pass resulting file paths to scientific commands. `biov exec
python` selects the scientific Python environment and does not include BioV.

## Remote execution

Set `execution_host` to an OpenSSH host alias and optionally set `ssh_config` and
`execution_cwd`. BioV sends `biov exec` and native arguments through SSH. The
remote BioV reads its own configuration and manifests; no files or local
settings are forwarded. Coordinates, native arguments, `--cwd` and
`--no-install` retain the same behavior remotely. Run setup on that execution host
with a local configuration. Setup rejects `execution_host` to prevent accidental local
installation for a remote deployment.

## Distribution and package sources

Pixi uses [BSD-3-Clause](https://github.com/prefix-dev/pixi/blob/v0.81.0/LICENSE).
BioV retains its notice when installing the manager. Distributors must also
satisfy applicable third-party binary and scientific package licenses.

The checked-in manifest names only conda-forge and Bioconda. Setup and execution
of declared environments use `--no-config` to ignore system/user Pixi
configuration; project-local Pixi
configuration is still respected. Custom manifests, locks and project settings
are operator-controlled: review their sources before distribution. Do not add
`defaults` or `repo.anaconda.com` to the shipped environment. Manager licensing,
repository access terms and individual package licenses are separate matters;
see [Anaconda's terms](https://www.anaconda.com/legal).

Native references: [pyproject.toml](https://pixi.prefix.dev/latest/python/pyproject_toml/),
[lock files](https://pixi.prefix.dev/latest/workspace/lockfile/),
[configuration](https://pixi.prefix.dev/latest/reference/pixi_configuration/),
[execution](https://pixi.prefix.dev/latest/reference/cli/pixi/run/).

## Planned tool lifecycle

This section is the broader design and acceptance target, not a claim that
every lifecycle operation is implemented. The upstream install/list/uninstall
adapters are described above; upgrades, general cache cleanup and broader tool
coverage remain gaps. See SPEC D6–D8.

The useful distinction from [uv's tool model](https://docs.astral.sh/uv/concepts/tools/)
is durable installation versus running without a persistent installation.
The latter can still reuse a disposable cached environment. uv also separates
[persistent tools, caches and command locations](https://docs.astral.sh/uv/reference/storage/).
BioV should expose those ownership and retention distinctions across its supported
biology tools, while continuing to use their native package managers and commands.

### First scope

Use a small set of upstream-managed tools on Linux x86_64, keeping the bundled
locked Pixi workflow separate. Add enough lifecycle control to inspect and reuse
environments, update explicitly, uninstall selected installations and clean
disposable content safely. Native upstream receipts/manifests are the installation
source of truth; content-addressed scientific workspaces have their own
provisioning/lock contract and are not the user-tool inventory.
Do not require every tool to support prebuilt binaries, containers and source
builds. Choose one verified route per tool and report its platform prerequisites.
Broader platform coverage and distributed management can follow demonstrated need.

Keep catalog discovery (what a tool does, its inputs and outputs) distinct from
installation discovery (what is actually installed here). Neither a successful
installation nor a zero exit status establishes scientific validity. Existing
identifier, native-file, coordinate, analysis-record and full-result reuse
contracts continue to apply.

### Broader acceptance targets

1. On a fresh supported host, prepare and run a selected supported tool without
   manually creating its environment. Report the actual source, version and
   environment; missing prerequisites and unsupported platforms are explicit.
2. Run the same locked tool from two unrelated analysis directories. Reuse the
   same ready environment without rebuilding; preserve each caller's directory,
   literal arguments and native exit status.
3. Distinguish durable installation from cached on-demand execution in inspection.
   Show owner, location, platform, resolved version/lock and readiness. Inspecting
   inventory does not install tools; a project-owned environment is not silently
   adopted as BioV-owned.
4. An exact pin selects that version or fails. Different resolved requirements
   stay isolated. Specify and test the precedence of an installed version, a
   cached version and a fresh resolution for unpinned requests. Refresh and
   upgrade are explicit; an upgrade respects constraints and does not rewrite a
   project lock silently.
5. Interrupted preparation never becomes a ready cache hit. A failed update leaves
   the previous working installation usable and exposes the original diagnostic.
6. Preview a selected uninstall or cache cleanup with ownership, paths and scope.
   After authorization, remove only eligible BioV-owned content; keep active runs,
   other retained installations, project environments, biological references,
   model weights and saved analysis inputs/outputs intact. A cleared on-demand
   environment can be recreated on the next run.
7. Storage locations are inspectable and user-configurable. No silent global PATH
   or shell-profile change, unrelated executable overwrite, or wrapper recursion
   occurs. Tool selection must not depend on remembering a build directory.
8. Re-run the existing representative two-step analysis and execution acceptance
   cases. Lifecycle changes preserve native data, scientific checks, complete
   result reuse and current command meanings.

The implemented local command surface is `install`, `list` and `uninstall`.
The remaining acceptance targets must be validated before extending this scope.
No tool registry service or MCP deployment-management surface is introduced.

## Standalone GOATOOLS enrichment

The bundled locked `goatools` environment pins GOATOOLS 1.6.5 and statsmodels
0.14.6, with its own Python 3.12 runtime. The separate uv tool installation pins
those same two packages but uses the paired BioV Python and upstream resolution
for the rest of its dependencies. GOATOOLS 1.6.5 imports
`multipletests` from statsmodels' old sandbox location; a real enrichment run
failed with statsmodels 0.15.0, so the compatible version is explicitly pinned
here. This does not alter the generic `python` or `txgnn` environments.
The example below exercises the existing Python-backed environment route.
The native install/direct-command route is described above; neither is an MCP operation. The Pixi `goatools` task delegates directly to the
[upstream 1.6.5 console entry point](https://github.com/tanghaibao/goatools/blob/v1.6.5/setup.cfg).
GO parsing and statistics remain in GOATOOLS and its scientific dependencies.

These are distinct facts:

- **Declared:** the shipped manifest and lock include `goatools`; this alone
  does not download packages or establish that the host can run it
- **Provisioned:** `biov python setup goatools` provisions the locked runtime;
  `biov install goatools` separately creates a uv tool environment and exposes
  its upstream console scripts; it does not provision the bundled locked workspace
- **Used:** `biov exec --no-install goatools ...` runs its upstream CLI, and a
  successful enrichment with independently checked outputs establishes use

For a clean isolated installation, set both roots. `BIOV_HOME` controls data
caches; the managed environment location is separately configured:

```bash
export BIOV_HOME="$PWD/.biov-go-data"
export BIOV_ENVIRONMENT_ROOT="$PWD/.biov-go-environments"
biov python setup goatools
biov exec --no-install goatools find_enrichment --help
```

With explicit, matching gene identifiers, population, annotations, and a local
ontology supplied by the caller:

```bash
biov exec --no-install goatools find_enrichment study.txt population.txt annotations.id2gos \
  --annofmt=id2gos --obo=ontology.obo --ns=BP --alpha=0.05 \
  --method=bonferroni,fdr_bh --pval=1 --pvalcalc=fisher_scipy_stats --outfile=results.tsv
```

`--pval` is an output filter applied after multiple-testing correction;
`--pval=1` retains all tested results, including non-significant terms. Default
GOATOOLS `is_a` count propagation is enabled. Choose these settings for the
scientific question rather than copying defaults blindly. BioV does not infer
identifier namespaces, choose annotations, select the background, or update
ontologies. Inputs and full native TSV results remain ordinary files.

From this repository, the default deterministic acceptance case invokes
`goatools` directly on PATH. First create the explicit installation and apply
the PATH instructions printed by install, if needed. Python compatibility setup
alone does not export this command:

```bash
biov install goatools
# Ensure environment_root/uv-tools/bin is on PATH in this shell.
goatools find_enrichment --help
python scripts/validate_goatools.py ./goatools-validation
```

The Python-backed manual analysis command above remains available as
`biov exec --no-install goatools ...`. For the acceptance script through native
provisioned execution instead of the exported command, pass
`--native-binary "$(command -v biov)"`; this uses `biov tools exec --no-install`
and requires its native provisioning receipt.

It uses only `tests/fixtures/goatools` inputs, runs the upstream GOATOOLS CLI, retains its
full TSV/stdout/stderr and settings with input/output SHA-256 values, and checks
the term set, study/population count ratios, and Fisher, Bonferroni and BH
values by independent exact integer arithmetic. The synthetic case tests three
BP hypotheses, including the propagated root, before filtering output.
`validation.json` verifies these selected outputs; it must be accompanied by
installation and version evidence to establish runtime identity. It is not
self-contained proof of which runtime executed the analysis.
Earlier acceptance of a clean cloud Linux-64 manager bootstrap, bundled locked
install, and this case passed with GOATOOLS 1.6.5, SciPy 1.18.1 and statsmodels 0.14.6.
This bounded fit case does not establish large-GO-database performance,
annotation/cache policy, support for other platforms, or biological validity.

The lock was regenerated by pinned Pixi 0.81.0 using its supported targeted
workflow, preserving every pre-existing environment's resolved package set:

```bash
biov pixi update --manifest-path src/biov/assets/environments/pyproject.toml \
  --environment goatools --no-config
```

## Native Rust locked-tool migration

The first native local setup-to-execution path supports the bundled `samtools`
and `goatools` environments on Linux x86_64. It uses the same unchanged native
Pixi manifest and lock as the Python route. All BioV control logic is Rust;
GOATOOLS itself still uses its locked Python runtime. The upstream programs keep
their native arguments and scientific behavior.

Build/install the standalone binary from this checkout:

```bash
cargo install --locked --path crates/biov-cli
biov tools --help
```

An existing Pixi 0.81.0 is required. The native route checks its reported version
and does not bootstrap/download a manager. Select it with `--pixi FILE` or
`BIOV_PIXI_BIN`; otherwise the matching managed location under the selected root
is preferred, then `pixi` on PATH. A configured incompatible executable fails
without switching sources. Existing `biov python setup goatools` can provision the
manager first if you use the Python route.

```bash
export BIOV_ENVIRONMENT_ROOT="$PWD/.biov-tool-environments"
export BIOV_PIXI_BIN="/path/to/existing/pixi"
biov tools inspect goatools
biov tools exec goatools find_enrichment --help
biov tools exec --no-install goatools find_enrichment --help
biov tools exec samtools --version
```

The executable example path is a placeholder for your matching local manager.
Without an explicit environment root, the native Linux default is
`$XDG_DATA_HOME/biov/environments`, or `$HOME/.local/share/biov/environments`.
`--environment-root` overrides `BIOV_ENVIRONMENT_ROOT`. Roots and manager/native
cache locations are separate; Pixi's native cache environment variables remain
available. `BIOV_HOME` is inherited without rewriting provider caches.
Native paths are literal OS paths: quoted `~` is not expanded. Use `$HOME` or
unquoted shell expansion. This intentionally differs from the Python route,
which applies `Path.expanduser()` to configured paths.

Default `tools exec` performs locked setup when there is no matching successful
native setup receipt and usable selected executable. A missing executable is
repaired from the unchanged lock: normal install is attempted first, then
`pixi reinstall --locked` if it is still missing. Redirected prefixes and unsafe
entrypoint symlinks fail before repair. `--no-install` never provisions: it
requires that receipt and a consistent Pixi prefix marker. Existing
Python-provisioned prefixes need a native `tools exec` (without `--no-install`) to establish the new receipt.
The separate `biov install` adapter does not create that receipt. Inspection works
without running the manager or creating a workspace and returns one JSON record on
stdout. Its `setup_recorded` status identifies recorded installation facts,
not independently verified package integrity, availability of the manager at
inspection time, biological validity or successful execution. Installation diagnostics use stderr; execution inherits native stdio and exit status
(including the conventional 128+signal exit code for signal termination).
Argument/configuration, setup, activation and pre-launch errors return 2 with
BioV diagnostics; an upstream program may also legitimately return 2. Exit status
alone does not distinguish these cases. The native status is preserved once the
selected executable has started.
During native `tools exec` CLI execution, direct SIGINT,
SIGTERM and SIGHUP to BioV are forwarded. BioV waits for its native child and
retains the workspace lock until that child exits. Noninteractive execution
forwards to a separate child process group, including ordinary descendants.
Interactive execution keeps the terminal's foreground group and forwards directly
to the selected child; terminal-generated signals retain normal group delivery.
Detached descendants, programs that ignore signals and SIGKILL are outside this
cooperative cleanup guarantee. This does not establish exactly-once delivery
across child startup or simultaneous user/group interrupts.

BioV options precede the tool name. An optional initial `--` after the name is
removed; everything else, including option-looking arguments, quoted strings,
empty strings and shell metacharacters, is passed as literal UTF-8 argv.
Non-UTF-8 arguments fail clearly before execution. The initial native route
reads typed `pixi shell-hook --as-is --json` activation, verifies its prefix,
project root/manifest, environment name and PATH prefix, then launches the exact
installed executable with OS argv, so a missing entry point cannot fall back to a same-name host executable.
No `pixi run` command-string parsing is used; zero native arguments stay zero,
including roots with spaces/apostrophes. Pixi already evaluates activation
scripts when computing its JSON environment; their exported variables are
preserved without sourcing the scripts a second time. Some punctuation-heavy roots remain
limited by pinned Pixi activation; failures or path expansions fail closed. It
intentionally bypasses the trivial bundled GOATOOLS shell task; native task
interpolation is not needed for these two entry points. `--cwd DIR`
selects the existing local analysis directory without moving the installation:

```bash
biov tools exec --no-install --cwd ./analysis goatools find_enrichment --help
```

This explicitly local interface does not read Python TOML configuration,
select project manifests, route to SSH, alter shell profiles,
or silently adopt external environments. `biov install` separately delegates
user-tool installation and command exposure to Pixi global/uv tool in dedicated
roots; its installations do not use this bundled workflow lock. Other bundled
tools, especially those requiring preparation tasks, still use the Python route.
Upgrades and general cache cleanup remain unimplemented. Workspace locks
coordinate native locked setup/runs only; Python/direct Pixi operations and hostile concurrent writers are outside that
coordination contract. A failed initial setup has no success receipt; repairing
a damaged prefix is not a transactional upgrade.

The fake-manager Rust acceptance tests verify pins, failure and busy-state
handling, literal argv/stdio/status, unchanged bundled bytes and cross-task reuse.
The existing deterministic GOATOOLS enrichment fixture can also be exercised
through this native entry point, with the same independent term/count/Fisher/
Bonferroni/BH checks. Neither test establishes general package/platform coverage.

For the independent GOATOOLS scientific acceptance case after native setup:

```bash
python scripts/validate_goatools.py ./native-go-validation \
  --native-binary "$(command -v biov)"
```

Earlier manual evidence for the separate bundled locked workflow from
2026-10-04 UTC on a Linux x86_64 cloud test executor
(Pixi 0.81.0, Rust 1.89.0): isolated native setup and this full enrichment check passed with
GOATOOLS 1.6.5, SciPy 1.18.1 and statsmodels 0.14.6. The test included a
quote-bearing output path. Default native setup/exec also installed and reported
Samtools 1.24 with HTSlib 1.24 from the unchanged bundled lock. These observations
cover those two selected distributions and this synthetic scientific example.

Reproduce real-manager regression coverage, distinct from fake-manager tests:

```bash
export BIOV_TEST_REAL_PIXI="/path/to/existing/pixi"
cargo test --locked -p biov-tools --lib real_pixi_tests -- --ignored --nocapture
```

The activation/argv test uses synthetic prefixes and sentinel executables; it
checks zero arguments under spaces/apostrophes, literal UTF-8 argv, native exit
status/signals and fail-closed handling of expanded activation paths. It does
not prove scientific package installation. The separate repair test installs
locked Samtools into its own temporary prefix, removes its executable, repairs
it and verifies unchanged lock bytes plus a successful native `--version` run.
Use native Pixi cache settings for offline/cached package availability. Both
opt-in tests passed with the pinned real manager on the same dated Linux host.
The GOATOOLS fixture above separately checks scientific enrichment results.

The new tools crate is workspace-coupled to authoritative bundled assets and
explicitly disables registry publishing. Supported native source builds use the
complete checkout or sdist; an installed binary runs outside the checkout.
