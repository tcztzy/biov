"""Acceptance tests for whole-script local and LSF execution."""

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

import biov.cli as cli
import biov.execution as execution_module
from biov.execution import (
    LSF_SUBMISSION_TIMEOUT_SECONDS,
    ExecutionSubmissionError,
    LsfSubmission,
    MissingExecutorError,
    execute_script,
)


def test_local_executor_runs_complete_script_with_arguments(tmp_path: Path) -> None:
    """Run the same ordinary Python script body in the local executor."""
    script = tmp_path / "analysis.py"
    output = tmp_path / "result.txt"
    script.write_text(
        "import pathlib, sys\npathlib.Path(sys.argv[1]).write_text(sys.argv[2])\n",
        encoding="utf-8",
    )

    result = execute_script(
        script,
        arguments=(str(output), "computed"),
        executor="local",
        python_executable=sys.executable,
        cwd=tmp_path,
    )

    assert isinstance(result, subprocess.CompletedProcess)
    assert result.returncode == 0
    assert output.read_text(encoding="utf-8") == "computed"
    assert result.args == (
        sys.executable,
        str(script.resolve()),
        str(output),
        "computed",
    )


def test_local_cli_propagates_script_exit_code(tmp_path: Path) -> None:
    """Make deterministic script failure visible as the CLI process status."""
    script = tmp_path / "failure.py"
    script.write_text("raise SystemExit(7)\n", encoding="utf-8")

    result = CliRunner().invoke(
        cli.app,
        ["run", "--executor", "local", str(script)],
    )

    assert result.exit_code == 7


def test_cli_forwards_arguments_after_separator(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Keep script options opaque to BioV when separated with ``--``."""
    script = tmp_path / "analysis.py"
    script.write_text("print('ok')\n", encoding="utf-8")
    seen: list[tuple[str, ...]] = []

    def fake_execute(*args, **kwargs):
        seen.append(kwargs["arguments"])
        return subprocess.CompletedProcess((sys.executable, str(script)), 0)

    monkeypatch.setattr(cli, "execute_script", fake_execute)

    result = CliRunner().invoke(
        cli.app,
        ["run", str(script), "--", "--input", "value with spaces"],
    )

    assert result.exit_code == 0
    assert seen == [("--input", "value with spaces")]


def test_lsf_executor_submits_argv_without_shell_and_returns_receipt(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Model accepted bsub output as a job receipt rather than completion."""
    script = tmp_path / "analysis.py"
    script.write_text("print('ok')\n", encoding="utf-8")
    calls: list[tuple[list[str], dict[str, object]]] = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="Job <4321> is submitted to queue <short>.\n",
            stderr="",
        )

    monkeypatch.setattr(execution_module.shutil, "which", lambda command: "/bin/bsub")
    monkeypatch.setattr(execution_module.subprocess, "run", fake_run)

    result = execute_script(
        script,
        arguments=("--input", "value with spaces"),
        executor="lsf",
        python_executable="/shared/venv/bin/python",
        cwd=tmp_path,
        queue="short",
        job_name="gc-analysis",
        stdout=tmp_path / "job.%J.out",
        stderr=tmp_path / "job.%J.err",
    )

    assert isinstance(result, LsfSubmission)
    assert result.job_id == 4321
    assert result.completed is False
    command, kwargs = calls[0]
    assert command == [
        "/bin/bsub",
        "-cwd",
        str(tmp_path.resolve()),
        "-q",
        "short",
        "-J",
        "gc-analysis",
        "-o",
        str(tmp_path / "job.%J.out"),
        "-e",
        str(tmp_path / "job.%J.err"),
        "/shared/venv/bin/python",
        str(script.resolve()),
        "--input",
        "value with spaces",
    ]
    assert kwargs["cwd"] == tmp_path.resolve()
    assert kwargs["capture_output"] is True
    assert kwargs["text"] is True
    assert kwargs["timeout"] == LSF_SUBMISSION_TIMEOUT_SECONDS
    assert "shell" not in kwargs or kwargs["shell"] is False


def test_lsf_executor_distinguishes_missing_rejected_and_unknown_receipt(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Expose stable errors for three different submission failure modes."""
    script = tmp_path / "analysis.py"
    script.write_text("print('ok')\n", encoding="utf-8")
    monkeypatch.setattr(execution_module.shutil, "which", lambda command: None)

    with pytest.raises(MissingExecutorError, match="bsub"):
        execute_script(script, executor="lsf")

    monkeypatch.setattr(execution_module.shutil, "which", lambda command: "/bin/bsub")
    monkeypatch.setattr(
        execution_module.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(
            returncode=2,
            stdout="",
            stderr="Queue is closed",
        ),
    )
    with pytest.raises(ExecutionSubmissionError, match="Queue is closed"):
        execute_script(script, executor="lsf")

    monkeypatch.setattr(
        execution_module.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(
            returncode=0,
            stdout="Submitted successfully",
            stderr="",
        ),
    )
    with pytest.raises(ExecutionSubmissionError, match="job ID"):
        execute_script(script, executor="lsf")


def test_lsf_executor_timeout_is_a_submission_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Bound a hung bsub instead of blocking the submitter forever."""
    script = tmp_path / "analysis.py"
    script.write_text("print('ok')\n", encoding="utf-8")
    monkeypatch.setattr(execution_module.shutil, "which", lambda command: "/bin/bsub")

    def hanging_run(command, **kwargs):
        raise subprocess.TimeoutExpired(command, kwargs["timeout"])

    monkeypatch.setattr(execution_module.subprocess, "run", hanging_run)

    with pytest.raises(ExecutionSubmissionError, match="did not respond"):
        execute_script(script, executor="lsf")


def test_lsf_cli_says_submitted_not_completed(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Keep the asynchronous LSF lifecycle explicit in user-visible output."""
    script = tmp_path / "analysis.py"
    script.write_text("print('ok')\n", encoding="utf-8")
    submission = LsfSubmission(
        job_id=4321,
        command=("bsub", "python", str(script)),
        message="Job <4321> is submitted to queue <normal>.",
    )
    monkeypatch.setattr(cli, "execute_script", lambda *args, **kwargs: submission)

    result = CliRunner().invoke(
        cli.app,
        ["run", "--executor", "lsf", str(script)],
    )

    assert result.exit_code == 0
    assert "Submitted LSF job 4321" in result.stdout
    assert "not completed" in result.stdout
    assert "completed successfully" not in result.stdout
