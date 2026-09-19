"""Whole-script execution in local and LSF environments."""

import os
import re
import shutil
import subprocess  # noqa: S404 - argv-only execution is this module's purpose
import sys
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

LSF_SUBMISSION_TIMEOUT_SECONDS = 60.0

_LSF_JOB_RECEIPT = re.compile(r"Job <(?P<job_id>[0-9]+)> is submitted\b")


class ExecutorKind(StrEnum):
    """Supported whole-script execution environments."""

    LOCAL = "local"
    LSF = "lsf"


class ExecutionError(RuntimeError):
    """Base class for expected whole-script execution failures."""


class MissingExecutorError(ExecutionError):
    """The selected local command or scheduler client is unavailable."""


class ExecutionSubmissionError(ExecutionError):
    """LSF rejected a job or returned no verifiable submission receipt."""


@dataclass(frozen=True, slots=True)
class LsfSubmission:
    """Accepted LSF job receipt; it is deliberately not a completion result."""

    job_id: int
    command: tuple[str, ...]
    message: str
    completed: bool = field(default=False, init=False)


def execute_script(
    script: Path,
    *,
    arguments: tuple[str, ...] = (),
    executor: ExecutorKind | str = ExecutorKind.LOCAL,
    python_executable: str | None = None,
    cwd: Path | None = None,
    queue: str | None = None,
    job_name: str | None = None,
    stdout: Path | None = None,
    stderr: Path | None = None,
) -> subprocess.CompletedProcess[bytes] | LsfSubmission:
    """Run or submit one complete ordinary Python analysis script.

    BioV does not inspect or rewrite the script. Calls to ``path`` occur
    inside the selected environment, so storage paths never cross the executor
    boundary.

    Args:
        script: Python script visible from the selected environment.
        arguments: Exact script argv following its path.
        executor: ``local`` or ``lsf``.
        python_executable: Interpreter command visible from the executor.
        cwd: Working directory inherited locally or pinned with LSF ``-cwd``.
        queue: Optional LSF queue.
        job_name: Optional LSF job name.
        stdout: Optional LSF stdout path.
        stderr: Optional LSF stderr path.

    Returns:
        Standard subprocess result or an asynchronous LSF submission receipt.

    Raises:
        ExecutionError: If validation, process launch, or submission fails.
        MissingExecutorError: If Python or the selected scheduler is unavailable.
        ExecutionSubmissionError: If LSF launch, acceptance, or receipt parsing fails.
    """
    try:
        selected = ExecutorKind(executor)
    except ValueError as error:
        raise ExecutionError(f"Unknown executor {executor!r}") from error
    try:
        script_path = Path(script).expanduser().resolve(strict=True)
    except OSError as error:
        raise ExecutionError(f"Script does not exist: {script}") from error
    if not script_path.is_file():
        raise ExecutionError(f"Script is not a regular file: {script}")
    candidate = Path.cwd() if cwd is None else cwd
    try:
        working_directory = candidate.expanduser().resolve(strict=True)
    except OSError as error:
        raise ExecutionError(
            f"Working directory does not exist: {candidate}"
        ) from error
    if not working_directory.is_dir():
        raise ExecutionError(f"Working directory is not a directory: {candidate}")
    if python_executable is None:
        python_executable = (
            os.getenv("BIOV_LSF_PYTHON") if selected is ExecutorKind.LSF else None
        ) or sys.executable
    if not python_executable or "\x00" in python_executable:
        raise MissingExecutorError("Python executable is empty or invalid")
    python_command = (python_executable, str(script_path), *arguments)
    if selected is ExecutorKind.LOCAL:
        try:
            return subprocess.run(  # noqa: S603
                python_command, cwd=working_directory, check=False
            )
        except FileNotFoundError as error:
            raise MissingExecutorError(
                f"Python executable is unavailable: {python_command[0]}"
            ) from error
        except OSError as error:
            raise ExecutionError(
                f"Could not start Python executable: {python_command[0]}"
            ) from error
    bsub = shutil.which("bsub")
    if bsub is None:
        raise MissingExecutorError("LSF executor requires bsub on PATH")
    command = [bsub, "-cwd", str(working_directory)]
    for flag, value, label in (("-q", queue, "queue"), ("-J", job_name, "job name")):
        if value is not None:
            if not value or any(character in value for character in "\x00\n\r"):
                raise ExecutionError(f"Invalid LSF {label}")
            command.extend((flag, value))
    if stdout is not None:
        command.extend(("-o", str(stdout)))
    if stderr is not None:
        command.extend(("-e", str(stderr)))
    command.extend(python_command)
    try:
        result = subprocess.run(  # noqa: S603
            command,
            cwd=working_directory,
            capture_output=True,
            text=True,
            check=False,
            timeout=LSF_SUBMISSION_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as error:
        raise ExecutionSubmissionError("LSF bsub did not respond in time") from error
    except FileNotFoundError as error:  # pragma: no cover - guarded by which
        raise MissingExecutorError("LSF bsub disappeared before submission") from error
    except OSError as error:
        raise ExecutionSubmissionError("Could not start LSF bsub") from error
    output = result.stdout.strip()
    error_output = result.stderr.strip()
    if result.returncode != 0:
        detail = (error_output or output or "no scheduler message")[:1000]
        raise ExecutionSubmissionError(f"LSF rejected job submission: {detail}")
    receipt = _LSF_JOB_RECEIPT.search(output)
    if receipt is None:
        detail = (output or error_output or "empty scheduler response")[:1000]
        raise ExecutionSubmissionError(
            f"LSF accepted no parseable numeric job ID: {detail}"
        )
    return LsfSubmission(
        job_id=int(receipt.group("job_id")),
        command=tuple(command),
        message=output,
    )


__all__ = [
    "LSF_SUBMISSION_TIMEOUT_SECONDS",
    "ExecutionError",
    "ExecutionSubmissionError",
    "ExecutorKind",
    "LsfSubmission",
    "MissingExecutorError",
    "execute_script",
]
