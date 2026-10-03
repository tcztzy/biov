"""Installed native FASTA-window MCP acceptance and offline Arrow portability.

Set BIOV_TEST_BINARY and BIOV_STANDALONE_PYTHON to run the synthetic gate.
The reader must contain only pyarrow==25.0.1, with BioV absent. An existing
GCF_000005845.2 package can additionally be selected with
BIOV_TEST_REFSEQ_PACKAGE; this test never downloads or changes original data.
"""

import asyncio
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
from collections import Counter
from contextlib import asynccontextmanager
from fractions import Fraction
from pathlib import Path

import pytest
from test_native_bundle_portability import _configured_executable
from test_prepared_fasta_portability import (
    _ACCESSION,
    _FASTA_PATH,
    _REFERENCE,
    _ROOT,
    _SEQUENCES,
    _environment,
    _hashes,
    _native_metadata,
    _offline_command,
    _register,
    _request,
    _synthetic_source,
)

_TIMEOUT = 120
_COLUMNS = [
    "sequence_id",
    "start",
    "end",
    "length",
    "is_full_window",
    "canonical_base_count",
    "gc_base_count",
    "gc_fraction",
    "weighted_gc_fraction",
]
_DTYPES = ["str", "i64", "i64", "i64", "bool", "i64", "i64", "f64", "f64"]
_IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def _metrics(sequence: str) -> dict[str, int | float | None]:
    """Return independent canonical and IUPAC-base-set reference metrics."""
    bases = sequence.upper()
    counts = Counter(bases)
    canonical = sum(counts[base] for base in "ACGT")
    gc = counts["G"] + counts["C"]
    weighted = sum(
        (
            Fraction(sum(base in "GC" for base in _IUPAC[code]), len(_IUPAC[code]))
            * count
            for code, count in counts.items()
        ),
        start=Fraction(),
    )
    return {
        "length": len(bases),
        "canonical_base_count": canonical,
        "gc_base_count": gc,
        "gc_fraction": gc / canonical if canonical else None,
        "weighted_gc_fraction": float(weighted / len(bases)) if bases else None,
    }


def _expected_rows(sequence_id: str, sequence: str, width: int) -> list[dict]:
    """Return all independently calculated zero-based half-open windows."""
    return [
        {
            "sequence_id": sequence_id,
            "start": start,
            "end": min(start + width, len(sequence)),
            "length": min(width, len(sequence) - start),
            "is_full_window": start + width <= len(sequence),
            **_metrics(sequence[start : start + width]),
        }
        for start in range(0, len(sequence), width)
    ]


def _assert_rows(actual: list[dict], expected: list[dict]) -> None:
    """Check complete typed records, with a narrow floating-point tolerance."""
    assert len(actual) == len(expected)
    for observed, reference in zip(actual, expected, strict=True):
        assert observed.keys() == reference.keys()
        for key, value in reference.items():
            if isinstance(value, float):
                assert math.isclose(observed[key], value, rel_tol=1e-14, abs_tol=1e-14)
            else:
                assert observed[key] == value, (key, observed, reference)


@asynccontextmanager
async def _mcp(binary: Path, data: Path, output: Path, store: Path | None):
    """Yield a real native MCP request function in a kernel-offline subprocess."""
    arguments = [
        str(binary),
        "mcp",
        "--data-root",
        str(data),
        "--output-root",
        str(output),
    ]
    if store is not None:
        arguments += ["--store-root", str(store)]
    process = await asyncio.create_subprocess_exec(
        *_offline_command(arguments),
        cwd=data.parent,
        env=_environment(),
        stdin=asyncio.subprocess.PIPE,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        limit=128 * 1024,
    )
    counter = 0

    async def request(method: str, params: dict) -> dict:
        nonlocal counter
        counter += 1
        assert process.stdin is not None and process.stdout is not None
        process.stdin.write(
            json.dumps(
                {
                    "jsonrpc": "2.0",
                    "id": counter,
                    "method": method,
                    "params": params,
                }
            ).encode()
            + b"\n"
        )
        await process.stdin.drain()
        while True:
            line = await asyncio.wait_for(process.stdout.readline(), _TIMEOUT)
            assert line, f"MCP exited before replying to {method}"
            reply = json.loads(line)
            assert reply["jsonrpc"] == "2.0", reply
            if "id" not in reply:
                assert "method" in reply, reply
                continue
            assert reply["id"] == counter and "error" not in reply, reply
            return reply["result"]

    async def call(name: str, arguments: dict, *, success: bool = True) -> dict:
        result = await asyncio.wait_for(
            request("tools/call", {"name": name, "arguments": arguments}),
            _TIMEOUT,
        )
        assert bool(result.get("isError")) is not success, (name, result)
        if success:
            assert isinstance(result.get("structuredContent"), dict), result
            return result["structuredContent"]
        assert result.get("content"), result
        return result

    try:
        initialized = await request(
            "initialize",
            {
                "protocolVersion": "2025-06-18",
                "capabilities": {},
                "clientInfo": {"name": "fasta-metrics-portability", "version": "1"},
            },
        )
        assert initialized["serverInfo"]["name"] == "biov-rs"
        assert process.stdin is not None
        process.stdin.write(b'{"jsonrpc":"2.0","method":"notifications/initialized"}\n')
        await process.stdin.drain()
        listed = await request("tools/list", {})
        tools = {tool["name"]: tool for tool in listed["tools"]}
        assert "dataset_fasta_windows" in tools
        schema = tools["dataset_fasta_windows"]["inputSchema"]
        assert set(schema["required"]) >= {
            "reference",
            "snapshot_id",
            "source_path",
            "recipe_id",
            "sequence_id",
            "window_size",
        }
        assert schema["properties"]["preview_rows"]["default"] == 5
        yield call
        process.stdin.close()
        stdout, stderr = await asyncio.wait_for(process.communicate(), _TIMEOUT)
        assert process.returncode == 0 and not stdout and not stderr, (stdout, stderr)
    finally:
        if process.returncode is None:
            process.kill()
            await asyncio.wait_for(process.communicate(), _TIMEOUT)


def _check_windows(result: dict, sequence_id: str, sequence: str, width: int) -> None:
    """Check a whole-sequence summary and the bounded first-window preview."""
    expected = _expected_rows(sequence_id, sequence, width)
    assert result["row_count"] == len(expected)
    assert result["schema"] == [
        {"name": name, "dtype": dtype}
        for name, dtype in zip(_COLUMNS, _DTYPES, strict=True)
    ]
    summary = result["sequence_summary"]
    assert summary["sequence_id"] == sequence_id
    for key, value in _metrics(sequence).items():
        if isinstance(value, float):
            assert math.isclose(summary[key], value, rel_tol=1e-14, abs_tol=1e-14)
        else:
            assert summary[key] == value, (key, summary, value)
    preview = result["preview"]
    count = min(5, len(expected))
    assert preview["returned_rows"] == count
    _assert_rows(
        [dict(zip(_COLUMNS, row, strict=True)) for row in preview["rows"]],
        expected[:count],
    )


def _copy_export(exported: dict, output: Path, destination: Path) -> None:
    """Copy exactly the four independently readable bundle files."""
    destination.mkdir()
    for key in ("execution_host_path", "record_path", "manifest_path", "readme_path"):
        path = Path(exported[key])
        assert path.is_file() and path.parent == output
        shutil.copyfile(path, destination / path.name)
    assert len(list(destination.iterdir())) == 4


async def _export_windows(
    binary: Path,
    data: Path,
    output: Path,
    store: Path,
    snapshot: dict,
    source_path: str,
    sequences: dict[str, str],
    width: int,
) -> tuple[dict, dict, str]:
    """Create full and filtered artifacts through actual MCP sequence/data tools.

    Returns:
        The full export, filtered export and historical filtered dataset handle.
    """
    async with _mcp(binary, data, output, store) as call:
        prepared = await call("prepared_fasta", _request(snapshot, source_path))
        base = {**_request(snapshot, source_path), "recipe_id": prepared["recipe_id"]}
        first = None
        for sequence_id, sequence in sequences.items():
            result = await call(
                "dataset_fasta_windows",
                {
                    **base,
                    "sequence_id": sequence_id,
                    "window_size": width,
                },
            )
            _check_windows(result, sequence_id, sequence, width)
            expected_rows = _expected_rows(sequence_id, sequence, width)
            if len(expected_rows) <= 50:
                complete = await call(
                    "dataset_preview",
                    {
                        "dataset_id": result["dataset_id"],
                        "preview_rows": 50,
                    },
                )
                _assert_rows(
                    [
                        dict(zip(_COLUMNS, row, strict=True))
                        for row in complete["preview"]["rows"]
                    ],
                    expected_rows,
                )
            if first is None:
                first = result
            else:
                await call("dataset_release", {"dataset_id": result["dataset_id"]})
        assert first is not None
        sequence_id = next(iter(sequences))
        for update in (
            {"recipe_id": "sha256-" + "0" * 64},
            {"sequence_id": "missing.0"},
            {"sequence_id": sequence_id.swapcase()},
            {"window_size": 0},
            {"window_size": 1024 * 1024 + 1},
            {"preview_rows": 51},
        ):
            await call(
                "dataset_fasta_windows",
                {
                    **base,
                    "sequence_id": sequence_id,
                    "window_size": width,
                    **update,
                },
                success=False,
            )
        prepared_path = store / "prepared" / prepared["recipe_id"]
        hidden_prepared = data.parent / "missing-preparation"
        prepared_path.rename(hidden_prepared)
        try:
            await call(
                "dataset_fasta_windows",
                {
                    **base,
                    "sequence_id": sequence_id,
                    "window_size": width,
                },
                success=False,
            )
            assert not prepared_path.exists(), (
                "Reading must not create a missing preparation"
            )
        finally:
            hidden_prepared.rename(prepared_path)
        native_fasta = store / snapshot["snapshot_path"] / "source" / source_path
        original_bytes = native_fasta.read_bytes()
        # Change one base in the disposable registered copy, never the import
        # or user's provider package. Failure must preserve corruption evidence.
        header_end = original_bytes.index(b"\n") + 1
        damaged_bytes = (
            original_bytes[:header_end]
            + (b"T" if original_bytes[header_end : header_end + 1] != b"T" else b"A")
            + original_bytes[header_end + 1 :]
        )
        native_fasta.write_bytes(damaged_bytes)
        try:
            damaged_store = _hashes(store)
            await call(
                "dataset_fasta_windows",
                {
                    **base,
                    "sequence_id": sequence_id,
                    "window_size": width,
                },
                success=False,
            )
            assert _hashes(store) == damaged_store, (
                "Corruption must not trigger silent repair"
            )
        finally:
            native_fasta.write_bytes(original_bytes)
        full = await call("dataset_export", {"dataset_id": first["dataset_id"]})
        origin = full["record"]["provenance"]["sequence_origin"]
        assert origin["reference"] == _REFERENCE
        assert origin["snapshot_id"] == snapshot["snapshot_id"]
        assert origin["recipe_id"] == prepared["recipe_id"]
        for field, path in (
            ("fai_sha256", "sequences.fai"),
            ("dictionary_sha256", "sequences.tsv"),
        ):
            assert (
                origin[field]
                == hashlib.sha256((prepared_path / path).read_bytes()).hexdigest()
            )
        assert (
            full["record"]["provenance"]["source_sha256"]
            == hashlib.sha256(
                (
                    store / snapshot["snapshot_path"] / "source" / source_path
                ).read_bytes()
            ).hexdigest()
        )
        selected = await call(
            "dataset_query",
            {
                "dataset_id": first["dataset_id"],
                "filter": {"column": "is_full_window", "op": "eq", "value": True},
                "sort": {"column": "start", "descending": True},
                "select": _COLUMNS,
                "preview_rows": 1,
            },
        )
        expected = [
            row
            for row in _expected_rows(sequence_id, sequences[sequence_id], width)
            if row["is_full_window"]
        ][::-1]
        assert selected["dataset_id"] != first["dataset_id"]
        assert selected["row_count"] == len(expected)
        assert selected["preview"]["returned_rows"] == min(1, len(expected))
        _assert_rows(
            [
                dict(zip(_COLUMNS, row, strict=True))
                for row in selected["preview"]["rows"]
            ],
            expected[:1],
        )
        subset = await call("dataset_export", {"dataset_id": selected["dataset_id"]})
        return full, subset, selected["dataset_id"]


async def _reopen_windows(
    binary: Path, output: Path, stale_id: str, subset: dict
) -> dict:
    """Reopen a persisted sequence artifact in a new process without raw input.

    Returns:
        The new export of a queried reopened artifact.
    """
    reopened_output = output.parent / "reopened-output"
    reopened_output.mkdir()
    async with _mcp(binary, output, reopened_output, None) as call:
        origin = subset["record"]["provenance"]["sequence_origin"]
        await call(
            "dataset_fasta_windows",
            {
                "reference": origin["reference"],
                "snapshot_id": origin["snapshot_id"],
                "source_path": subset["record"]["provenance"]["source"],
                "recipe_id": origin["recipe_id"],
                "sequence_id": origin["sequence_id"],
                "window_size": origin["window_size"],
            },
            success=False,
        )
        await call("dataset_preview", {"dataset_id": stale_id}, success=False)
        reopened = await call(
            "dataset_reopen",
            {
                "record_path": Path(subset["record_path"]).name,
                "preview_rows": 1,
            },
        )
        assert reopened["dataset_id"] != stale_id
        assert reopened["row_count"] == subset["record"]["row_count"]
        assert reopened["schema"] == subset["record"]["schema"]
        assert (
            reopened["source_sha256"] == subset["record"]["provenance"]["source_sha256"]
        )
        assert (
            reopened["scientific_metadata"] == subset["record"]["scientific_metadata"]
        )
        queried = await call(
            "dataset_query",
            {
                "dataset_id": reopened["dataset_id"],
                "sort": {"column": "start", "descending": True},
                "preview_rows": 1,
            },
        )
        assert queried["row_count"] == reopened["row_count"]
        result = await call("dataset_export", {"dataset_id": queried["dataset_id"]})
        assert result["record"]["reopen_verification"] is not None
        assert (
            result["record"]["provenance"]["sequence_origin"]
            == subset["record"]["provenance"]["sequence_origin"]
        )
        return result


def _binary_and_reader() -> tuple[Path, Path]:
    """Return explicit installed executables, skipping only an unconfigured gate."""
    if not any(
        name in os.environ for name in ("BIOV_TEST_BINARY", "BIOV_STANDALONE_PYTHON")
    ):
        pytest.skip("Set BIOV_TEST_BINARY and BIOV_STANDALONE_PYTHON for acceptance")
    return (
        _configured_executable("BIOV_TEST_BINARY"),
        _configured_executable("BIOV_STANDALONE_PYTHON"),
    )


def _move_and_read(
    tmp_path: Path,
    reader_python: Path,
    source: Path,
    store: Path,
    output: Path,
    full: dict,
    subset: dict,
    reopened: dict,
    mode: str,
) -> None:
    """Remove historical dependencies and analyze complete moved IPC with PyArrow."""
    with tempfile.TemporaryDirectory(prefix="independent-fasta-metrics-") as directory:
        consumer = Path(directory).resolve()
        assert not consumer.is_relative_to(_ROOT) and not consumer.is_relative_to(
            tmp_path
        )
        _copy_export(full, output, consumer / "full")
        _copy_export(subset, output, consumer / "selected")
        reopened_output = Path(reopened["record_path"]).parent
        _copy_export(reopened, reopened_output, consumer / "reopened")
        reader = consumer / "read_fasta_metrics.py"
        shutil.copyfile(_ROOT / "tests/_fasta_metrics_reader.py", reader)
        shutil.rmtree(source)
        shutil.rmtree(store)
        shutil.rmtree(output)
        shutil.rmtree(reopened_output)
        assert not any(
            path.exists() for path in (source, store, output, reopened_output)
        )
        before = _hashes(consumer)
        completed = subprocess.run(
            _offline_command([str(reader_python), "-I", str(reader), mode]),
            cwd=consumer,
            env=_environment(),
            capture_output=True,
            text=True,
            timeout=_TIMEOUT,
            check=False,
        )
        assert completed.returncode == 0, (completed.stdout, completed.stderr)
        assert "standalone FASTA metrics filter and summary passed" in completed.stdout
        assert _hashes(consumer) == before


@pytest.mark.parametrize("newline", [b"\n", b"\r\n"], ids=["lf", "crlf"])
def test_fasta_metrics_synthetic_mcp_arrow_move(tmp_path: Path, newline: bytes) -> None:
    """Check mixed-case/IUPAC metrics and complete native table reuse without BioV."""
    binary, reader = _binary_and_reader()
    data, store, output = (tmp_path / name for name in ("imports", "store", "output"))
    source = data / "synthetic"
    fasta = _synthetic_source(source, newline)
    extra_id, extra_sequence = "iupac-All:15", "acgtRYSWKMBDHVNryswkmbdhvn"
    fasta.write_bytes(
        fasta.read_bytes()
        + newline
        + f">{extra_id} complete IUPAC fixture".encode()
        + newline
        + extra_sequence.encode()
    )
    _native_metadata(source, fasta)
    sequences = {**_SEQUENCES, extra_id: extra_sequence}
    store.mkdir()
    output.mkdir()
    source_before = _hashes(source)
    snapshot = _register(binary, tmp_path, store, source)
    full, subset, stale_id = asyncio.run(
        _export_windows(
            binary,
            data,
            output,
            store,
            snapshot,
            _FASTA_PATH,
            sequences,
            4,
        )
    )
    assert _hashes(source) == source_before
    store_before = _hashes(store)
    reopened = asyncio.run(_reopen_windows(binary, output, stale_id, subset))
    assert _hashes(store) == store_before and _hashes(source) == source_before
    _move_and_read(
        tmp_path, reader, source, store, output, full, subset, reopened, "synthetic"
    )


def test_fasta_metrics_real_refseq_mcp_arrow_move(tmp_path: Path) -> None:
    """Check explicit 10 kb windows of an existing E. coli package without downloads."""
    if "BIOV_TEST_REFSEQ_PACKAGE" not in os.environ:
        pytest.skip("Set BIOV_TEST_REFSEQ_PACKAGE for the existing RefSeq package")
    binary, reader = _binary_and_reader()
    original = Path(os.environ["BIOV_TEST_REFSEQ_PACKAGE"]).absolute()
    assert original.is_dir()
    before = _hashes(original)
    data, store, output = (tmp_path / name for name in ("imports", "store", "output"))
    source = data / "real-refseq"
    shutil.copytree(original, source)
    source_path = f"ncbi_dataset/data/{_ACCESSION}/{_ACCESSION}_ASM584v2_genomic.fna"
    lines = (source / source_path).read_text(encoding="ascii").splitlines()
    assert lines[0].split()[0] == ">NC_000913.3"
    assert not any(line.startswith(">") for line in lines[1:])
    sequence = "".join(lines[1:])
    assert len(sequence) == 4_641_652
    store.mkdir()
    output.mkdir()
    try:
        snapshot = _register(binary, tmp_path, store, source)
        full, subset, stale_id = asyncio.run(
            _export_windows(
                binary,
                data,
                output,
                store,
                snapshot,
                source_path,
                {"NC_000913.3": sequence},
                10_000,
            )
        )
        assert full["record"]["row_count"] == 465
        assert subset["record"]["row_count"] == 464
        reopened = asyncio.run(_reopen_windows(binary, output, stale_id, subset))
        _move_and_read(
            tmp_path, reader, source, store, output, full, subset, reopened, "real"
        )
    finally:
        assert _hashes(original) == before, "Original downloaded package changed"


async def _reject_row_budget(
    binary: Path, data: Path, output: Path, store: Path, snapshot: dict
) -> None:
    """Prove preflight rejections preserve disk and every dataset handle slot."""
    async with _mcp(binary, data, output, store) as call:
        prepared = await call("prepared_fasta", _request(snapshot))
        base = {
            **_request(snapshot),
            "recipe_id": prepared["recipe_id"],
            "sequence_id": "budget.1",
        }
        before = _hashes(store)
        for _ in range(2):
            rejected = await call(
                "dataset_fasta_windows",
                {
                    **base,
                    "window_size": 1,
                    "preview_rows": 0,
                },
                success=False,
            )
            assert "budget" in json.dumps(rejected).lower()
        assert _hashes(store) == before
        assert _hashes(output) == {}
        successful = await call(
            "dataset_fasta_windows",
            {
                **base,
                "window_size": 10_000,
                "preview_rows": 0,
            },
        )
        assert successful["row_count"] == 41
        assert successful["preview"]["returned_rows"] == 0
        # There are sixteen session dataset slots. Failed preflights must not
        # consume one: this successful sequence result leaves exactly fifteen.
        for _ in range(15):
            opened = await call("dataset_open", {"path": "tiny.csv", "preview_rows": 0})
            assert opened["row_count"] == 1
        await call("dataset_open", {"path": "tiny.csv"}, success=False)
        assert _hashes(store) == before and _hashes(output) == {}


def test_fasta_metrics_row_budget_rejection_is_atomic(tmp_path: Path) -> None:
    """Reject a 400,001-row request before allocating a frame or publishing a handle."""
    binary, _ = _binary_and_reader()
    data, store, output = (tmp_path / name for name in ("imports", "store", "output"))
    source = data / "budget"
    fasta = _synthetic_source(source)
    sequence = "ACGT" * 100_000 + "A"
    fasta.write_text(
        ">budget.1 generated row-budget fixture\n"
        + "\n".join(
            sequence[start : start + 80] for start in range(0, len(sequence), 80)
        ),
        encoding="ascii",
    )
    _native_metadata(source, fasta)
    (data / "tiny.csv").write_text("id\n001\n", encoding="ascii")
    store.mkdir()
    output.mkdir()
    before = _hashes(source)
    snapshot = _register(binary, tmp_path, store, source)
    asyncio.run(_reject_row_budget(binary, data, output, store, snapshot))
    assert _hashes(source) == before
