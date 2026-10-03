"""Run package coordinates on the configured execution host."""

import json
import os
import shutil
import subprocess  # noqa: S404 - explicit argv is the execution interface
from pathlib import Path

from .config import settings
from .environments import (
    declared_environments,
    environment_manifest,
    pixi_command,
    provision_environment,
)
from .remote import run_remote


def _has_task(manager: tuple[str, ...], manifest: Path, name: str) -> bool | None:
    """Ask Pixi whether this environment resolves a same-name task.

    Returns:
        Whether a workspace or feature task provides the entry point, or None
        when Pixi cannot report its tasks, so the answer stays unknown rather
        than being guessed.
    """
    try:
        result = subprocess.run(  # noqa: S603
            [
                *manager,
                "task",
                "--no-config",
                "list",
                "--json",
                "--manifest-path",
                str(manifest),
            ],
            capture_output=True,
            check=True,
        )
        listing = json.loads(result.stdout)
    except (OSError, KeyError, TypeError, ValueError, subprocess.SubprocessError):
        # Introspection is advisory: it must never replace the entry's own status.
        return None
    return any(
        task["name"] == name
        for record in listing
        if record["environment"] == name
        for group in (record, *record["features"])
        for task in group["tasks"]
    )


def _missing_entry_reason(
    manager: tuple[str, ...],
    manifest: Path,
    name: str,
    *,
    working_directory: Path,
    environment: dict[str, str],
    known_task: bool | None,
) -> str | None:
    """Explain a 127 exit that no declared entry point accounts for.

    Returns:
        Error message when neither a same-name task nor a same-name executable
        exists, otherwise None, including when Pixi cannot be asked at all.
    """
    try:
        if known_task is None:
            known_task = _has_task(manager, manifest, name)
        if known_task is not False:
            return None
        activation = subprocess.run(  # noqa: S603
            [
                *manager,
                "shell-hook",
                "--no-config",
                "--as-is",
                "--json",
                "--manifest-path",
                str(manifest),
                "--environment",
                name,
            ],
            cwd=working_directory,
            env=environment,
            capture_output=True,
            check=True,
        )
        path = json.loads(activation.stdout)["environment_variables"]["PATH"]
        if shutil.which(name, path=path) is not None:
            return None
        features = declared_environments()[name]
    except (OSError, KeyError, TypeError, ValueError, subprocess.SubprocessError):
        return None
    instruction = (
        f"[tool.pixi.feature.{json.dumps(features[0])}.tasks] "
        f'{json.dumps(name)} = "<command>"'
        if features
        else f"select a feature defining a {name!r} task"
    )
    return (
        f"declared environment {name!r} declares no {name!r} task or "
        f"executable; add one in {manifest}: {instruction}; "
        f"if the environment is not installed, run `biov setup {name}`"
    )


def run_software(
    tool: str,
    arguments: tuple[str, ...] = (),
    *,
    cwd: Path | None = None,
    install: bool = True,
) -> subprocess.CompletedProcess[bytes]:
    """Run a package coordinate locally or through SSH.

    Bare names prefer declared locked Pixi environments, then temporary conda
    environments. The conda, pypi and npm prefixes select Pixi, uv and npx respectively;
    conda names also reuse a declared environment when one exists.

    Args:
        tool: NAME or SOURCE:NAME package coordinate.
        arguments: Native arguments passed unchanged.
        cwd: Working directory on the execution host.
        install: Provision a declared Pixi environment first; False skips
            installation and preparation and refuses temporary package environments.

    Returns:
        Process status with inherited stdin, stdout and stderr.

    Raises:
        ValueError: If a package name is invalid or temporary installation is disabled.
        NotADirectoryError: If the working directory is not a directory.
        FileNotFoundError: If a required manager or declared entry point is absent.
    """
    cwd = settings.execution_cwd if cwd is None else cwd
    if settings.execution_host is not None:
        command = ["biov", "exec"]
        if not install:
            command.append("--no-install")
        if cwd is not None:
            command.extend(("--cwd", str(cwd)))
        return run_remote(
            settings.execution_host,
            (*command, "--", tool, *arguments),
            ssh_config=settings.ssh_config,
        )
    source, separator, name = tool.partition(":")
    if not separator or source not in {"conda", "pypi", "npm"}:
        source, name = "conda", tool
    if not name or name.startswith("-"):
        raise ValueError(
            f"{source}: requires a package name that does not start with '-'"
        )
    working_directory = (
        (Path.cwd() if cwd is None else cwd).expanduser().resolve(strict=True)
    )
    if not working_directory.is_dir():
        raise NotADirectoryError(working_directory)
    command = [name, *arguments]
    declared = source == "conda" and name in declared_environments()
    workspace: tuple[tuple[str, ...], Path] | None = None
    task: bool | None = None
    if not declared and not install:
        raise ValueError(
            f"--no-install cannot run temporary package {tool!r}; omit it or "
            "use a declared Pixi environment"
        )
    if source == "conda":
        # Resolve the manager before the command's working directory applies, and
        # provision with that same executable instead of resolving it twice.
        manager = pixi_command()
        if not Path(manager[0]).is_file():
            raise FileNotFoundError(
                f"No Pixi executable at {manager[0]}; run `biov setup` to install "
                "the pinned manager"
            )
        if declared:
            manifest = (
                provision_environment(name, executable=manager[0])
                if install
                else environment_manifest()
            )
            if not manifest.is_file():
                raise FileNotFoundError(
                    f"declared environment {name!r} has no initialized workspace "
                    f"at {manifest}; run `biov setup {name}` before --no-install"
                )
            manifest = manifest.resolve(strict=True)
            workspace = manager, manifest
            if any("'" in argument for argument in arguments):
                task = _has_task(manager, manifest, name)
                if task:
                    # Pixi 0.81 quotes task arguments without escaping single
                    # quotes. Escape those quotes inside its surrounding quotes.
                    command = [
                        name,
                        *(arg.replace("'", "'\"'\"'") for arg in arguments),
                    ]
            command = [
                *manager,
                "run",
                "--no-config",
                # On-demand provisioning installs from the lock, so the run itself
                # installs from that lock too; --no-install skips both
                # installation and lock updates.
                "--frozen" if install else "--as-is",
                "--manifest-path",
                str(manifest),
                "--environment",
                name,
                "--",
                *command,
            ]
        else:
            command = [*manager, "exec", "-s", name, "--", *command]
    else:
        executable, *options = (
            ("uv", "tool", "run") if source == "pypi" else ("npx", "--yes")
        )
        installed = shutil.which(executable)
        if installed is None:
            installation = (
                "https://docs.astral.sh/uv/getting-started/installation/"
                if source == "pypi"
                else "https://docs.npmjs.com/downloading-and-installing-node-js-and-npm"
            )
            raise FileNotFoundError(
                f"{tool} requires {executable} on PATH; install it: {installation}"
            )
        command = [str(Path(installed).resolve()), *options, *command]
    cache = settings.home.expanduser().resolve()
    environment = {
        **os.environ,
        "BIOV_HOME": str(cache),
        "BIOV_CACHE_HTTP": str(settings.cache_http).lower(),
    }
    result = subprocess.run(  # noqa: S603
        command,
        cwd=working_directory,
        env=environment,
        check=False,
    )
    if workspace is not None and result.returncode == 127:
        manager, manifest = workspace
        reason = _missing_entry_reason(
            manager,
            manifest,
            name,
            working_directory=working_directory,
            environment=environment,
            known_task=task,
        )
        if reason is not None:
            raise FileNotFoundError(reason)
    return result
