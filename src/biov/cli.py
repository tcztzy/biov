"""BioV command-line interface."""

import json
import subprocess  # noqa: S404 - pass native arguments to the resolved manager
from pathlib import Path
from typing import Annotated

import typer
from typer.core import TyperCommand

from .config import ConfigFileError, select_config_file
from .environments import pixi_command, setup_environment
from .execution import (
    ExecutionError,
    ExecutorKind,
    LsfSubmission,
    execute_script,
)
from .registry import RegistryAssetError, RegistryUpdateResult, update_registry_asset
from .software import run_software

app = typer.Typer(no_args_is_help=True)


@app.callback()
def main(
    config: Annotated[
        Path | None,
        typer.Option(
            "--config",
            help="Application TOML file; overrides BIOV_CONFIG for this invocation.",
        ),
    ] = None,
) -> None:
    """BioV utilities.

    Raises:
        typer.Exit: With a stable status when the selected configuration is unusable.
    """
    if config is None:
        return
    try:
        select_config_file(config)
    except ConfigFileError as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(2) from error


@app.command("mcp")
def run_mcp() -> None:
    """Run the BioV MCP server over standard input/output."""
    from .mcp import main as run_mcp_server

    run_mcp_server()


@app.command("analyze")
def analyze(
    request: Annotated[Path, typer.Argument(help="Analysis request JSON file")],
) -> None:
    """Run one recorded analysis in a declared local Pixi environment.

    Raises:
        typer.Exit: With status 1 for a failed run or 2 for an unusable request.
    """
    from .analysis import AnalysisRequest, run_analysis

    try:
        result = run_analysis(AnalysisRequest.model_validate_json(request.read_text()))
    except (OSError, ValueError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(2) from error
    typer.echo(json.dumps(result, ensure_ascii=False))
    if result["status"] != "succeeded":
        raise typer.Exit(1)


@app.command("inspect-analysis")
def inspect_saved_analysis(
    record: Annotated[
        str, typer.Argument(help="Saved analysis record path or file URI")
    ],
) -> None:
    """Inspect a saved run without executing or resubmitting it.

    Raises:
        typer.Exit: With status 2 if the record or its results cannot be read.
    """
    from .analysis import inspect_analysis

    try:
        result = inspect_analysis(record)
    except (OSError, ValueError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(2) from error
    typer.echo(json.dumps(result, ensure_ascii=False))


@app.command("setup")
def setup_software(
    tool: Annotated[
        str | None,
        typer.Argument(help="Declared Pixi environment to install"),
    ] = None,
    archive: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            dir_okay=False,
            help="Pinned official Pixi archive for offline manager installation",
        ),
    ] = None,
    all_environments: Annotated[
        bool,
        typer.Option("--all", help="Install every declared scientific environment"),
    ] = False,
    update_lock: Annotated[
        bool,
        typer.Option(
            help="Regenerate the explicitly selected project's lock before installation"
        ),
    ] = False,
) -> None:
    """Install Pixi, optionally installing a locked scientific environment.

    Raises:
        typer.Exit: With a stable status when setup cannot proceed.
    """
    try:
        installed = setup_environment(
            tool,
            archive=archive,
            all_environments=all_environments,
            update_lock=update_lock,
        )
    except ValueError as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(2) from error
    typer.echo(installed)


class _NativePixiCommand(TyperCommand):
    """Keep Pixi argv intact instead of asking Typer to parse its options."""

    def parse_args(self, ctx, args: list[str]) -> list[str]:
        """Let Typer initialize the context without consuming native tokens.

        Returns:
            Arguments Typer itself still has to process.
        """
        remaining = super().parse_args(ctx, [])
        ctx.params["arguments"] = list(args)
        return remaining


@app.command("pixi", cls=_NativePixiCommand, add_help_option=False)
def run_pixi(
    arguments: Annotated[
        list[str] | None,
        typer.Argument(help="Native Pixi arguments; passed unchanged"),
    ] = None,
) -> None:
    """Run the resolved Pixi manager with native arguments.

    Raises:
        typer.Exit: With Pixi's exit status, or a stable launch error status.
    """
    command = pixi_command()
    if not Path(command[0]).is_file():
        typer.echo(
            f"No Pixi executable at {command[0]}; run `biov setup` to install the "
            "pinned manager",
            err=True,
        )
        raise typer.Exit(2)
    result = subprocess.run([*command, *(arguments or ())], check=False)  # noqa: S603
    raise typer.Exit(result.returncode)


@app.command(
    "exec",
    context_settings={"ignore_unknown_options": True, "allow_interspersed_args": False},
)
def execute_software(
    tool: Annotated[
        str, typer.Argument(help="[conda:|pypi:|npm:]NAME; default source: conda")
    ],
    arguments: Annotated[
        list[str] | None,
        typer.Argument(help="Native arguments; passed unchanged after NAME"),
    ] = None,
    cwd: Annotated[
        Path | None,
        typer.Option(help="Analysis working directory on the execution host"),
    ] = None,
    no_install: Annotated[
        bool,
        typer.Option(
            "--no-install",
            help=(
                "Run only what is already installed: skip installing the selected "
                "Pixi environment and its prepare task"
            ),
        ),
    ] = False,
) -> None:
    """Run a package coordinate; place BioV options before NAME.

    Declared Pixi environments use their lock and preparation task. Other
    names use temporary package environments, defaulting to conda.

    Raises:
        typer.Exit: With the native program's exit status or a launch error status.
    """
    native = tuple(arguments or ())
    if native[:1] == ("--",):
        native = native[1:]
    try:
        result = run_software(tool, native, cwd=cwd, install=not no_install)
    except (OSError, ValueError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(2) from error
    raise typer.Exit(result.returncode)


@app.command("run")
def run_script(
    script: Annotated[
        Path,
        typer.Argument(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            help="Ordinary Python analysis script to execute",
        ),
    ],
    arguments: Annotated[
        list[str] | None,
        typer.Argument(
            help="Arguments passed unchanged to the script; use -- before options"
        ),
    ] = None,
    executor: Annotated[
        ExecutorKind,
        typer.Option(help="Whole-script execution environment"),
    ] = ExecutorKind.LOCAL,
    python_executable: Annotated[
        str | None,
        typer.Option(
            "--python",
            help=(
                "Python executable visible in the selected environment; LSF also "
                "reads BIOV_LSF_PYTHON"
            ),
        ),
    ] = None,
    cwd: Annotated[
        Path | None,
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            help="Execution working directory; defaults to the current directory",
        ),
    ] = None,
    queue: Annotated[
        str | None,
        typer.Option(help="LSF queue"),
    ] = None,
    job_name: Annotated[
        str | None,
        typer.Option(help="LSF job name"),
    ] = None,
    stdout: Annotated[
        Path | None,
        typer.Option(help="LSF stdout path; %J is replaced by the job ID"),
    ] = None,
    stderr: Annotated[
        Path | None,
        typer.Option(help="LSF stderr path; %J is replaced by the job ID"),
    ] = None,
) -> None:
    """Run a complete Python script locally or submit it to LSF.

    Raises:
        typer.Exit: With the local script status or a stable launch error status.
    """
    try:
        result = execute_script(
            script,
            arguments=tuple(arguments or ()),
            executor=executor,
            python_executable=python_executable,
            cwd=cwd,
            queue=queue,
            job_name=job_name,
            stdout=stdout,
            stderr=stderr,
        )
    except ExecutionError as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(2) from error
    if not isinstance(result, LsfSubmission):
        raise typer.Exit(result.returncode)
    typer.echo(
        f"Submitted LSF job {result.job_id}; submission accepted, job not completed."
    )


@app.command("update")
def update_registry(
    output: Annotated[
        Path | None,
        typer.Option(
            help="Registry asset destination; defaults to the packaged asset."
        ),
    ] = None,
    force: Annotated[
        bool,
        typer.Option(
            help="Replace the asset even when its response bytes are unchanged."
        ),
    ] = False,
    timeout: Annotated[
        float,
        typer.Option(help="Network timeout in seconds", min=0.001),
    ] = 60,
) -> None:
    """Synchronize the complete identifiers.org registry asset.

    Raises:
        typer.Exit: With a stable status when validation or fetching fails.
    """
    try:
        result = update_registry_asset(output, force=force, timeout=timeout)
    except (RegistryAssetError, OSError) as error:
        typer.echo(str(error), err=True)
        raise typer.Exit(2) from error
    if result.status == "unchanged":
        typer.echo(
            f"identifiers.org registry is already current; skipped {result.output}"
        )
        return
    namespaces = result.asset["payload"]["namespaces"]
    typer.echo(f"Updated {len(namespaces)} namespaces in {result.output}")


__all__ = [
    "RegistryUpdateResult",
    "app",
    "main",
    "run_mcp",
    "run_script",
    "update_registry",
]
