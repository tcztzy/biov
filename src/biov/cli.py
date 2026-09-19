"""BioV command-line interface."""

from pathlib import Path
from typing import Annotated

import typer

from .execution import (
    ExecutionError,
    ExecutorKind,
    LsfSubmission,
    execute_script,
)
from .registry import RegistryAssetError, RegistryUpdateResult, update_registry_asset

app = typer.Typer(no_args_is_help=True)


@app.callback()
def main() -> None:
    """BioV utilities."""


@app.command("mcp")
def run_mcp() -> None:
    """Run the BioV MCP server over standard input/output."""
    from .mcp import main as run_mcp_server

    run_mcp_server()


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


@app.command("update-identifiers-registry")
def update_identifiers_registry(
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
    "update_identifiers_registry",
]
