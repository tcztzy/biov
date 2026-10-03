"""Persist and inspect ordinary Python analyses in declared Pixi environments."""

import hashlib
import json
import mimetypes
import os
import re
import shutil
import stat
import subprocess  # noqa: S404 - the caller explicitly supplies analysis code
import sys
from importlib.metadata import version
from itertools import islice
from pathlib import Path
from typing import Any, Literal
from urllib.parse import quote, unquote, urlsplit
from uuid import uuid4

import pandas as pd
from Bio import SeqIO
from pydantic import BaseModel, ConfigDict, Field, JsonValue, field_validator

from .config import settings
from .environments import manifest_source, pixi_command, resolve_environment

PREVIEW_ROWS = 2
MAX_RESPONSE_BYTES = 32768
MAX_RESOURCE_BYTES = 1048576
_ENGINE_FILES = {
    "record.json",
    "execution.json",
    "code.py",
    "inputs.json",
    "parameters.json",
    "checks.json",
    "stdout.log",
    "stderr.log",
    "inputs",
    "record.tmp",
    "execution.tmp",
    "inputs.tmp",
    "parameters.tmp",
}
_READABLE_FILES = {
    "record.json",
    "execution.json",
    "code.py",
    "checks.json",
    "stdout.log",
    "stderr.log",
}


class AnalysisRequest(BaseModel):
    """One Python script, its explicit file inputs and its declared outputs."""

    model_config = ConfigDict(extra="forbid")

    name: str = Field(min_length=1, max_length=120)
    environment: str = Field(min_length=1, max_length=80)
    code: str = Field(min_length=1, max_length=131072)
    inputs: dict[str, str] = Field(default_factory=dict, max_length=16)
    outputs: dict[str, Literal["csv", "fasta", "file"]] = Field(
        min_length=1, max_length=16
    )
    parameters: dict[str, JsonValue] = Field(default_factory=dict)
    requirements: list[str] = Field(default_factory=list)

    @field_validator("inputs", "outputs")
    @classmethod
    def safe_names(cls, value: dict[str, Any]) -> dict[str, Any]:
        """Reject ambiguous names, engine files and directory traversal.

        Returns:
            The validated mapping, without renaming caller files.

        Raises:
            ValueError: If a name is unsafe, reserved or duplicated ignoring case.
        """
        for name in value:
            if (
                re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,119}", name) is None
                or ".." in name
                or name.casefold() in _ENGINE_FILES
            ):
                raise ValueError(f"Unsafe or reserved analysis filename: {name!r}")
        if len({name.casefold() for name in value}) != len(value):
            raise ValueError("Analysis filenames must also be distinct ignoring case")
        return value


class _Check(BaseModel):
    """A check actually reported by the supplied analysis script."""

    model_config = ConfigDict(extra="forbid", strict=True)
    name: str = Field(min_length=1)
    passed: bool
    detail: str


def _root() -> Path:
    """Return the configured persistent result directory."""
    return settings.analysis_root.expanduser().resolve()


def _local_path(reference: str) -> Path:
    """Read an ordinary local path or local file URI without resolving symlinks.

    Returns:
        Absolute local path.

    Raises:
        ValueError: If a URI names another host or storage scheme.
    """
    if "://" in reference:
        uri = urlsplit(reference)
        if (
            uri.scheme != "file"
            or uri.netloc not in {"", "localhost"}
            or uri.query
            or uri.fragment
        ):
            raise ValueError("Managed analysis accepts local paths or local file URIs")
        reference = unquote(uri.path)
    return Path(os.path.abspath(Path(reference).expanduser()))


def _digest(path: Path) -> str:
    """Return the SHA-256 of the file's complete bytes."""
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def _save(path: Path, record: dict[str, Any]) -> None:
    """Atomically update a record with one writer in its private run directory."""
    temporary = path.with_suffix(".tmp")
    with temporary.open("w", encoding="utf-8") as output:
        json.dump(record, output, ensure_ascii=False, allow_nan=False, indent=2)
        output.flush()
        os.fsync(output.fileno())
    temporary.replace(path)


def _record_path(reference: str) -> Path:
    """Validate a record's location within the current execution context.

    Returns:
        Existing record path.

    Raises:
        ValueError: If the reference is outside this server's results root.
    """
    path = _local_path(reference)
    if path.name != "record.json" or path.parent.parent != _root() or path.is_symlink():
        raise ValueError("Record reference does not belong to this analysis root")
    if path.resolve(strict=True).parent.parent != _root():
        raise ValueError("Record reference escapes its analysis root")
    return path


def _output_file(directory: Path, name: str, metadata: dict[str, Any]) -> Path:
    """Require an unchanged, independent regular output file.

    Returns:
        Validated output path.

    Raises:
        ValueError: If the output is replaced, linked or modified.
    """
    AnalysisRequest.safe_names({name: "file"})
    path = directory / name
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ValueError(f"Output is not an independent regular file: {name}")
    if info.st_size != metadata["size"] or _digest(path) != metadata["sha256"]:
        raise ValueError(f"Saved output identity mismatch: {name}")
    return path


def _source(reference: str) -> tuple[Path, str | None]:
    """Require a regular input and verify prior analysis outputs before reuse.

    Returns:
        Input path and its required prior-output digest, when applicable.

    Raises:
        ValueError: If a previous result is incomplete or has changed.
        FileNotFoundError: If the input is missing or is not a regular file.
    """
    path = _local_path(reference)
    if path.is_relative_to(_root()) or path.resolve().is_relative_to(_root()):
        record = json.loads(_record_path(str(path.parent / "record.json")).read_text())
        if record["status"] != "succeeded" or path.name not in record["outputs"]:
            raise ValueError(
                "Only completed registered outputs can be reused as inputs"
            )
        metadata = record["outputs"][path.name]
        return _output_file(path.parent, path.name, metadata), metadata["sha256"]
    if not path.is_file():
        raise FileNotFoundError(f"Analysis input is not a regular file: {path}")
    return path, None


def _snapshot(source: Path, target: Path) -> dict[str, Any]:
    """Copy a stable source to the private read-only input used by the script.

    Returns:
        Source and snapshot identity metadata.

    Raises:
        ValueError: If the source changes during copying.
    """

    def identity(info: os.stat_result) -> tuple[int, ...]:
        return (
            info.st_dev,
            info.st_ino,
            info.st_size,
            info.st_mtime_ns,
            info.st_ctime_ns,
        )

    with source.open("rb") as original, target.open("xb") as copied:
        before = identity(os.fstat(original.fileno()))
        shutil.copyfileobj(original, copied)
        if before != identity(os.fstat(original.fileno())) or before != identity(
            source.stat()
        ):
            raise ValueError(f"Input changed while copying: {source}")
    target.chmod(0o400)
    return {
        "source": source.as_uri(),
        "snapshot": target.as_uri(),
        "size": target.stat().st_size,
        "sha256": _digest(target),
    }


def _preview(path: Path, kind: str) -> dict[str, Any]:
    """Read the first records using mature parsers and bound displayed content.

    Returns:
        Preview metadata with explicit omissions and unscanned totals.
    """
    if kind == "csv":
        frame = pd.read_csv(path, nrows=PREVIEW_ROWS + 1)
        visible = frame.iloc[:PREVIEW_ROWS, :8]
        rows = json.loads(visible.to_json(orient="values"))
        shortened = False
        for row in rows:
            for index, value in enumerate(row):
                if isinstance(value, str) and len(value) > 128:
                    row[index] = value[:128] + "…"
                    shortened = True
        return {
            "columns": [str(name)[:128] for name in visible.columns],
            "dtypes": [str(dtype) for dtype in visible.dtypes],
            "rows": rows,
            "total_rows": len(frame) if len(frame) <= PREVIEW_ROWS else None,
            "total_columns": len(frame.columns),
            "order": "file",
            "truncated": len(frame) > PREVIEW_ROWS
            or len(frame.columns) > 8
            or shortened
            or any(len(str(name)) > 128 for name in visible.columns),
        }
    if kind == "fasta":
        with path.open() as source:
            sequences = list(islice(SeqIO.parse(source, "fasta"), PREVIEW_ROWS + 1))
        return {
            "records": [
                {
                    "id": item.id[:128],
                    "description": item.description[:128],
                    "length": len(item),
                    "sequence_prefix": str(item.seq[:80]),
                }
                for item in sequences[:PREVIEW_ROWS]
            ],
            "total_records": len(sequences) if len(sequences) <= PREVIEW_ROWS else None,
            "order": "file",
            "truncated": len(sequences) > PREVIEW_ROWS
            or any(len(item) > 80 or len(item.description) > 128 for item in sequences),
        }
    return {"content_preview": None, "reason": "Generic file; use the complete result"}


def _download(path: Path) -> str | None:
    """Return an operator-configured existing storage URL, without uploading."""
    if settings.analysis_base_url is None:
        return None
    return (
        str(settings.analysis_base_url).rstrip("/")
        + "/"
        + quote(path.relative_to(_root()).as_posix())
    )


def _summary(record: dict[str, Any]) -> dict[str, Any]:
    """Bound the whole default response while retaining its complete record.

    Returns:
        JSON-compatible response no larger than the default response ceiling.
    """
    result = {
        key: record[key]
        for key in (
            "record",
            "name",
            "status",
            "stage",
            "exit_code",
            "inputs",
            "outputs",
            "checks",
            "diagnostic",
            "logs",
            "environment",
            "parameters",
            "requirements",
        )
    }
    result = json.loads(json.dumps(result, ensure_ascii=False))
    omitted: list[str] = []
    if len(result["diagnostic"]) > 2048:
        result["diagnostic"] = result["diagnostic"][:2048]
        omitted.append("diagnostic remainder")
        result["omitted"] = omitted

    def size() -> int:
        return len(json.dumps(result, ensure_ascii=False, allow_nan=False).encode())

    # Leave room for the MCP result envelope without duplicating summary text.
    ceiling = MAX_RESPONSE_BYTES - 512
    if size() > ceiling:
        result["omitted"] = omitted
        for name, output in reversed(list(result["outputs"].items())):
            output["preview"] = {"omitted": "response size limit"}
            omitted.append(f"outputs.{name}.preview")
            if size() <= ceiling:
                return result
        for field in (
            "checks",
            "requirements",
            "parameters",
            "inputs",
            "environment",
            "outputs",
            "logs",
        ):
            result[field] = (
                {"name": record["environment"]["name"]}
                if field == "environment"
                else []
                if field in {"checks", "requirements"}
                else {}
            )
            omitted.append(field)
            if size() <= ceiling:
                return result
    return result


def inspect_analysis(record: str) -> dict[str, Any]:
    """Inspect saved facts without resubmitting or finalizing an unconfirmed run.

    File-access errors and invalid references or content identities propagate.

    Args:
        record: File URI or local path of the saved run record.

    Returns:
        Bounded run summary; querying a failed run is a successful query.

    """
    path = _record_path(record)
    saved = json.loads(path.read_text())
    if saved["status"] == "succeeded":
        for name, metadata in saved["outputs"].items():
            output = _output_file(path.parent, name, metadata)
            metadata["download_url"] = _download(output)
    elif saved["status"] != "failed":
        saved["status"] = "unknown"
        saved["diagnostic"] = (
            "No saved terminal result; execution status is unconfirmed. No resubmission was attempted."
        )
        evidence = path.parent / "execution.json"
        if evidence.is_file():
            saved["environment"]["execution"] = json.loads(evidence.read_text())
    return _summary(saved)


def read_analysis_file(uri: str) -> tuple[bytes, str]:
    """Read one registered complete result or diagnostic within the resource cap.

    File-access errors propagate if a registered file is missing or inaccessible.

    Args:
        uri: Ordinary file URI in the configured analysis root.

    Returns:
        Complete bytes and their MIME type.

    Raises:
        ValueError: If the reference, identity or response size is invalid.
    """
    path = _local_path(uri)
    record = json.loads(_record_path(str(path.parent / "record.json")).read_text())
    metadata = None
    if path.name in record["outputs"] and record["status"] == "succeeded":
        metadata = record["outputs"][path.name]
        path = _output_file(path.parent, path.name, record["outputs"][path.name])
    elif path.name not in _READABLE_FILES:
        raise ValueError("File is not a completed registered result or diagnostic")
    if path.is_symlink() or not path.is_file():
        raise ValueError("Analysis resource is not a regular file")
    with path.open("rb") as source:
        content = source.read(MAX_RESOURCE_BYTES + 1)
    if len(content) > MAX_RESOURCE_BYTES:
        location = _download(path)
        raise ValueError(
            f"Complete file exceeds the {MAX_RESOURCE_BYTES}-byte resource limit; "
            + (
                f"download from {location}"
                if location
                else "configure a reachable BIOV_ANALYSIS_BASE_URL for complete retrieval"
            )
        )
    if metadata is not None and (
        len(content) != metadata["size"]
        or hashlib.sha256(content).hexdigest() != metadata["sha256"]
    ):
        raise ValueError(f"Saved output changed during retrieval: {path.name}")
    return content, mimetypes.guess_type(path.name)[0] or "application/octet-stream"


def run_analysis(request: AnalysisRequest) -> dict[str, Any]:
    """Run a declared Python analysis with persistent inputs, outputs and facts.

    File-access errors propagate if the persistent run record cannot be saved.

    Args:
        request: Validated script, parameters and explicit file contract.

    Returns:
        Bounded summary, including failed or unconfirmed execution states.

    Raises:
        ValueError: If remote managed execution is selected but unsupported.
    """
    if settings.execution_host is not None:
        raise ValueError("Managed analysis currently supports local execution only")
    directory = _root() / uuid4().hex
    directory.mkdir(mode=0o700, parents=True)
    record_path = directory / "record.json"
    record: dict[str, Any] = {
        "record": record_path.as_uri(),
        "name": request.name,
        "status": "launch_intent",
        "stage": None,
        "exit_code": None,
        "inputs": {},
        "outputs": {},
        "checks": [],
        "diagnostic": "",
        "logs": {
            name: (directory / f"{name}.log").as_uri() for name in ("stdout", "stderr")
        },
        "environment": {"name": request.environment},
        "cwd": str(directory),
        "requirements": request.requirements,
        "parameters": request.parameters,
        "declared_outputs": request.outputs,
    }
    _save(record_path, record)
    for name in ("stdout", "stderr"):
        (directory / f"{name}.log").touch()
    stage = "data_checks"
    launcher_started = False
    confirmed_completion = False
    try:
        (directory / "inputs").mkdir()
        for name, reference in request.inputs.items():
            source, expected = _source(reference)
            snapshot = _snapshot(source, directory / "inputs" / name)
            if expected is not None and snapshot["sha256"] != expected:
                raise ValueError(f"Saved output changed before snapshot: {source}")
            record["inputs"][name] = snapshot
        _save(
            directory / "inputs.json",
            {
                name: str(_local_path(item["snapshot"]))
                for name, item in record["inputs"].items()
            },
        )
        _save(directory / "parameters.json", request.parameters)
        (directory / "code.py").write_text(request.code, encoding="utf-8")
        record["code"] = {
            "uri": (directory / "code.py").as_uri(),
            "sha256": _digest(directory / "code.py"),
        }
        stage = "launch"
        resolve_environment(request.environment)
        manager = pixi_command()
        manifest = manifest_source()
        manager_version = subprocess.run(  # noqa: S603
            [*manager, "--version"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        record["environment"].update(
            manifest=str(manifest),
            manifest_sha256=_digest(manifest),
            lock_sha256=_digest(manifest.with_name("pixi.lock")),
            pixi_version=manager_version,
            biov_driver_version=version("biov"),
            driver_interpreter=sys.executable,
        )
        command = [
            sys.executable,
            "-c",
            "from biov.cli import app; app()",
            "exec",
            "--cwd",
            str(directory),
            request.environment,
            str(Path(__file__).with_name("_analysis_worker.py")),
            str(directory),
        ]
        record["command"] = command
        _save(record_path, record)
        with (
            (directory / "stdout.log").open("wb") as stdout,
            (directory / "stderr.log").open("wb") as stderr,
        ):
            launcher = subprocess.Popen(  # noqa: S603
                command,
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
            )
            launcher_started = True
            record["launcher_pid"] = launcher.pid
            _save(record_path, record)
            launcher_code = launcher.wait()
        record["launcher_exit_code"] = launcher_code
        evidence_path = directory / "execution.json"
        if not evidence_path.exists():
            record.update(
                status="unknown",
                stage="launch",
                diagnostic="Launcher exited without scientific-process confirmation; execution status is unknown.",
            )
        else:
            evidence = json.loads(evidence_path.read_text())
            record["environment"]["execution"] = evidence
            if evidence["status"] == "failed":
                record.update(
                    status="failed",
                    stage=evidence["stage"],
                    diagnostic=evidence["diagnostic"],
                )
            elif evidence["status"] != "completed":
                record.update(
                    status="unknown",
                    diagnostic="Scientific execution has no confirmed completion; no resubmission was attempted.",
                )
            else:
                confirmed_completion = True
                record["exit_code"] = evidence["exit_code"]
                stage = "data_checks"
                for name, metadata in record["inputs"].items():
                    if _digest(directory / "inputs" / name) != metadata["sha256"]:
                        raise ValueError(f"Execution input identity mismatch: {name}")
                checks = directory / "checks.json"
                if checks.exists():
                    raw_checks = json.loads(checks.read_text())
                    if not isinstance(raw_checks, list):
                        raise ValueError(
                            "checks.json must contain a list of actual checks"
                        )
                    record["checks"] = [
                        _Check.model_validate(item).model_dump() for item in raw_checks
                    ]
                failed = [
                    item["name"] for item in record["checks"] if not item["passed"]
                ]
                if failed:
                    raise ValueError("Failed declared checks: " + ", ".join(failed))
                if evidence["exit_code"]:
                    record.update(
                        status="failed",
                        stage="program",
                        diagnostic=f"Scientific program exited with status {evidence['exit_code']}; inspect the complete logs.",
                    )
                else:
                    outputs = {}
                    for name, kind in request.outputs.items():
                        path = directory / name
                        info = path.lstat()
                        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                            raise ValueError(
                                f"Output is not an independent regular file: {name}"
                            )
                        outputs[name] = {
                            "uri": path.as_uri(),
                            "download_url": _download(path),
                            "format": kind,
                            "size": info.st_size,
                            "sha256": _digest(path),
                            "preview": _preview(path, kind),
                        }
                    record.update(status="succeeded", outputs=outputs)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        record.update(
            status="unknown"
            if launcher_started and not confirmed_completion
            else "failed",
            stage=stage,
            diagnostic=str(error),
        )
    _save(record_path, record)
    return _summary(record)
