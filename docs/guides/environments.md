# Scientific environments

BioV runs in the Python installation that already hosts the application. Pixi
manages separate scientific environments; BioV does not install a second copy
of itself. Environment definitions use native `[tool.pixi.*]` tables in
`src/biov/assets/environments/pyproject.toml`, with resolved versions in the
adjacent `pixi.lock`. The root project and development tools use uv.

## New execution host

```bash
biov setup
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
biov setup --archive /path/to/pixi-aarch64-apple-darwin.tar.gz
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
biov setup samtools
biov setup --all
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
`biov setup` after Pixi installs dependencies. They are not replaced by similarly
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
BIOV_ENVIRONMENT_MANIFEST=/path/to/pyproject.toml biov setup --update-lock
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
biov setup python
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
biov setup vina-meeko
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

This section is a design and acceptance target, not a list of commands available
today. Current setup and execution behavior is documented above. The next step is
to make BioV own tool environments across tasks, so neither the researcher nor an
agent has to remember where each tool was installed or built. See SPEC D6–D8.

The useful distinction from [uv's tool model](https://docs.astral.sh/uv/concepts/tools/)
is durable installation versus running without a persistent installation.
The latter can still reuse a disposable cached environment. uv also separates
[persistent tools, caches and command locations](https://docs.astral.sh/uv/reference/storage/).
BioV should expose those ownership and retention distinctions across its supported
biology tools, while continuing to use their native package managers and commands.

### First scope

Use a small set of existing locked Pixi tools on Linux x86_64. Add enough lifecycle
control to find and inspect environments, reuse them across directories, update
explicitly, uninstall selected installations and clean disposable content safely.
Existing manifest/lock content-addressed workspaces are a starting point; they
are not, by themselves, a complete installed-tool inventory or cleanup policy.
Do not require every tool to support prebuilt binaries, containers and source
builds. Choose one verified route per tool and report its platform prerequisites.
Broader platform coverage and distributed management can follow demonstrated need.

Keep catalog discovery (what a tool does, its inputs and outputs) distinct from
installation discovery (what is actually installed here). Neither a successful
installation nor a zero exit status establishes scientific validity. Existing
identifier, native-file, coordinate, analysis-record and full-result reuse
contracts continue to apply.

### Acceptance cases before implementation

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

Finalize exact command names, status fields and test fixtures in T56 before
implementation. This specification does not introduce a new CLI alias, package
format, tool registry service or MCP deployment-management surface.
