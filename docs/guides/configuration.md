# Configuration

Set environment variables before importing BioV or starting its CLI/MCP server.
BioV does not automatically load a project `.env` file.

BioV also reads `config.toml` from `platformdirs.user_config_path("biov")`.
`BIOV_CONFIG` selects another application TOML file, and the CLI's root-level
`--config PATH` option overrides `BIOV_CONFIG` for one invocation:

```console
biov --config /shared/project/biov.toml exec samtools view -c sample.bam
```

Both selections must name an existing TOML file. A missing path, a directory, or
an empty value fails with one error naming `BIOV_CONFIG` and `--config`; the
platform default file may be absent. Environment variables override file
settings. Deployment choices are private application settings, not MCP resources
or tool parameters.

```toml
execution_host = "compute"
execution_cwd = "/shared/project"
```

The installed SSH client resolves `compute` through its standard SSH
configuration. Omit `execution_host` for local execution. In either case,
callers use the same `biov exec [SOURCE:]NAME ARGS...` interface. Manifest and cache
settings are read on the execution host; local settings are not forwarded.
A bare name defaults to conda and selects a declared Pixi environment first,
otherwise temporary Pixi execution.
`conda:NAME` uses a declared environment when present, otherwise temporary Pixi
execution; `pypi:NAME` and `npm:NAME` always use their native package managers.
Declared environments use their lock and preparation tasks; temporary sources
resolve packages without that lock. There are no software-alias or default-runtime
settings. Missing uv or npx fails with installation instructions and no fallback.

See the [execution guide](environments.md) for package sources, `--no-install`
and remote-path behavior.

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `BIOV_HOME` | Platform cache dir | Custom cache directory |
| `BIOV_CACHE_HTTP` | `True` | Enable HTTP caching |
| `BIOV_CONFIG` | User config directory / `config.toml` | Application TOML file; override for one invocation with `--config PATH` |
| `BIOV_EXECUTION_HOST` | Unset | Private SSH host alias; unset means local execution |
| `BIOV_SSH_CONFIG` | OpenSSH defaults | Optional native SSH configuration file |
| `BIOV_EXECUTION_CWD` | Execution host's working directory | Default analysis working directory |
| `BIOV_ENVIRONMENT_MANIFEST` | Bundled `pyproject.toml` and `pixi.lock` | Explicit project manifest path on the execution host |
| `BIOV_ENVIRONMENT_ROOT` | Platform user data directory / `biov/environments` | Managed Pixi and scientific environments; used only when no configured or `PATH` Pixi matches the pinned version |
| `BIOV_PIXI_BIN` | Unset | Explicit Pixi executable; otherwise a matching `pixi` on `PATH`, otherwise the managed copy |
| `BIOV_ANALYSIS_ROOT` | Platform user data directory / `biov/results` | Persistent managed-analysis outputs and records; separate from disposable caches |
| `BIOV_ANALYSIS_BASE_URL` | Unset | Existing HTTP(S) storage prefix serving the results root, for client downloads without a shared filesystem; no query, fragment or embedded credentials |

See [managed analysis](analysis.md) for the first acceptance case and verification
status. Setting `BIOV_ANALYSIS_BASE_URL` does not start a server or publish files.

`BIOV_MAX_FILE_BYTES` optionally limits each downloaded or decompressed provider
file in bytes. It is unset by default because scientific file sizes vary widely.
An exceeded limit fails before the file is published to the artifact cache.

`BIOV_HOME` selects the cache root; `~` is expanded and relative paths are
resolved when settings load. If unset, BioV uses
`$XDG_CACHE_HOME/biov`, then the platform's user cache directory. Choose a local
directory or mounted disk in your shell or MCP launch configuration; no machine's
mount path is stored in the repository.

Identifier files live under `<BIOV_HOME>/artifacts/`. Ordinary fsspec
`filecache` reads use `BIOV_HOME` as their default directory. An existing fsspec
configuration (including `FSSPEC_FILECACHE`) is preserved, and an explicit
`cache_storage` argument takes precedence.

```python
from biov.config import settings

print(settings.config)  # application TOML file this process reads
print(settings.cache_http)
```
