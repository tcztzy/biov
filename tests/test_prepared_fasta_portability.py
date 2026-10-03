"""Installed prepared-FASTA acceptance with independent, offline scientific readers.

Set BIOV_TEST_BINARY and BIOV_STANDALONE_FASTA for the default synthetic tests.
The latter must be an isolated environment containing only pysam==0.24.1,
biopython==1.88 and numpy==2.5.3. No BioV Python package is imported here or by
that reader. BIOV_TEST_REFSEQ_PACKAGE additionally enables the already downloaded
GCF_000005845.2 package; these tests never download data or change originals.

A portable export is the exact selected native snapshot plus its prepared
folder, preserving their store-relative layout. A prepared folder alone does
not contain the raw FASTA. The independent reader starts with that folder only,
uses its provenance/README, and receives no CLI response or historical path.
"""

import ctypes
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

import pytest

_ROOT = Path(__file__).resolve().parents[1]
_TIMEOUT = 120
_REFERENCE = "refseq.gcf:GCF_000005845.2"
_ACCESSION = "GCF_000005845.2"
_FASTA_PATH = f"ncbi_dataset/data/{_ACCESSION}/genome.fna"
_ENVIRONMENT = ("BIOV_TEST_BINARY", "BIOV_STANDALONE_FASTA")
_SEQUENCES = {
    "chrAlpha.1": "ACGTacgtnNRYacGTACGTnnnnACGTacGTACGTA",
    "scaffold|two:2": "ttGGCCaaNNacgtACG",
    "plasmid-C.7": "n",
    "001": "aCgT",
}


def _configured_path(name: str, *, executable: bool = False) -> Path:
    """Return an explicit installed executable or existing fixture directory."""
    value = os.environ.get(name)
    assert value, f"Installed acceptance configuration is missing: {name}"
    path = Path(value).absolute()
    if executable:
        assert path.is_file() and os.access(path, os.X_OK), (name, path)
    else:
        assert path.is_dir(), (name, path)
    return path


def _binary() -> Path:
    """Return the installed binary, skipping only an unconfigured checkout."""
    if not any(name in os.environ for name in _ENVIRONMENT):
        pytest.skip("Set BIOV_TEST_BINARY and BIOV_STANDALONE_FASTA for acceptance")
    return _configured_path("BIOV_TEST_BINARY", executable=True)


def _offline_command(arguments: list[str]) -> list[str]:
    """Return a command using the existing child-only kernel network sandbox."""
    assert sys.platform == "linux", "Offline acceptance needs Linux + libseccomp"
    ctypes.CDLL("libseccomp.so.2")
    return [
        sys.executable,
        "-I",
        str(_ROOT / "tests/_native_storage_sandbox.py"),
        *arguments,
    ]


def _environment() -> dict[str, str]:
    """Return an environment without shell paths or checkout imports."""
    return {
        "PATH": "",
        "PYTHONPATH": "",
        "PYTHONNOUSERSITE": "1",
        "POLARS_MAX_THREADS": "2",
    }


def _run_offline(
    arguments: list[str], cwd: Path, *, success: bool = True
) -> subprocess.CompletedProcess[str]:
    """Return a bounded child result with success or rejection checked."""
    result = subprocess.run(
        _offline_command(arguments),
        cwd=cwd,
        env=_environment(),
        capture_output=True,
        text=True,
        timeout=_TIMEOUT,
        check=False,
    )
    if success:
        assert result.returncode == 0, (arguments[:3], result.stdout, result.stderr)
    else:
        assert result.returncode != 0, "Corrupt prepared/source data was accepted"
        assert not result.stdout.strip(), "Failure must not emit a success result"
        assert result.stderr.strip(), "Rejected preparation needs a diagnostic"
    return result


def _cli(
    binary: Path,
    cwd: Path,
    store: Path,
    request: dict[str, Any],
    *,
    source: Path | None = None,
    success: bool = True,
) -> dict[str, Any] | None:
    """Return registration/preparation JSON, or None for an expected rejection."""
    with tempfile.TemporaryDirectory(prefix="fasta-request-", dir=cwd) as directory:
        request_path = Path(directory) / "request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        if source is None:
            arguments = [str(binary), "prepared", "fasta"]
        else:
            arguments = [str(binary), "storage", "register"]
        arguments += ["--store-root", str(store), "--request-file", str(request_path)]
        if source is not None:
            arguments += ["--source-root", str(source.parent)]
        result = _run_offline(arguments, cwd, success=success)
        return json.loads(result.stdout) if success else None


def _hashes(root: Path) -> dict[str, tuple[int, str]]:
    """Return streamed file identities without following source symlinks."""
    identities = {}
    for path in sorted(root.rglob("*")):
        assert not path.is_symlink(), path
        if path.is_file():
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            identities[path.relative_to(root).as_posix()] = (
                path.stat().st_size,
                digest,
            )
    return identities


def _relative_path(store: Path, name: str) -> Path:
    """Return a validated, directly readable confined store-relative path."""
    assert name and not Path(name).is_absolute() and "\\" not in name, name
    assert all(part not in {"", ".", ".."} for part in name.split("/")), name
    path = store / name
    assert path.exists() and path.resolve().is_relative_to(store.resolve()), name
    assert not path.is_symlink(), name
    return path


def _native_metadata(source: Path, fasta: Path) -> None:
    """Describe generated bytes using the native catalog/checksum conventions."""
    data = source / "ncbi_dataset/data"
    catalog = data / "dataset_catalog.json"
    catalog.write_text(
        json.dumps(
            {
                "apiVersion": "V2",
                "assemblies": [
                    {
                        "accession": _ACCESSION,
                        "files": [
                            {
                                "filePath": fasta.relative_to(data).as_posix(),
                                "fileType": "GENOMIC_NUCLEOTIDE_FASTA",
                                "uncompressedLengthBytes": str(fasta.stat().st_size),
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    checksums = []
    for path in (catalog, fasta):
        with path.open("rb") as stream:
            digest = hashlib.file_digest(
                stream, lambda: hashlib.md5(usedforsecurity=False)
            ).hexdigest()
        checksums.append(f"{digest}  {path.relative_to(source).as_posix()}\n")
    (source / "md5sum.txt").write_text("".join(checksums), encoding="utf-8")


def _synthetic_source(source: Path, newline: bytes = b"\n") -> Path:
    """Return a generated native FASTA with mixed case, wrapping and exact IDs."""
    fasta = source / _FASTA_PATH
    fasta.parent.mkdir(parents=True)
    (source / "README.md").write_text("Synthetic fixture; not provider data.\n")
    lines = []
    for name, sequence in _SEQUENCES.items():
        lines.append(f">{name} descriptive text is not the sequence ID".encode())
        lines.extend(
            sequence[start : start + 7].encode() for start in range(0, len(sequence), 7)
        )
    # No final newline tests the last-record byte boundary independently of LF/CRLF.
    fasta.write_bytes(newline.join(lines))
    _native_metadata(source, fasta)
    return fasta


def _register(binary: Path, cwd: Path, store: Path, source: Path) -> dict[str, Any]:
    """Return the registered snapshot of a complete disposable native package."""
    result = _cli(
        binary,
        cwd,
        store,
        {
            "source_path": source.name,
            "requested_ref": _REFERENCE,
            "canonical_ref": _REFERENCE,
            "declaration": {"provider": "refseq"},
        },
        source=source,
    )
    assert result is not None
    return result["snapshot"]


def _request(
    snapshot: dict[str, Any], source_path: str = _FASTA_PATH
) -> dict[str, str]:
    """Return a request selecting an exact reference, snapshot and native file."""
    return {
        "reference": _REFERENCE,
        "snapshot_id": snapshot["snapshot_id"],
        "source_path": source_path,
    }


def _check_result(
    result: dict[str, Any], store: Path, snapshot: dict[str, Any]
) -> Path:
    """Return the verified prepared folder after checking ordinary output paths."""
    assert result["reference"] == _REFERENCE
    assert result["snapshot_id"] == snapshot["snapshot_id"]
    assert result["recipe_id"].startswith("sha256-")
    assert len(result["recipe_id"]) == 71
    paths = {}
    for field in ("fasta", "fai", "dictionary", "provenance", "readme"):
        item = result[field]
        path = _relative_path(store, item["relative_path"])
        assert path.is_file() and Path(item["execution_host_path"]) == path
        paths[field] = path
    prepared = store / "prepared" / result["recipe_id"]
    assert paths["fai"] == prepared / "sequences.fai"
    assert paths["dictionary"] == prepared / "sequences.tsv"
    assert paths["provenance"] == prepared / "provenance.json"
    assert paths["readme"] == prepared / "README.md"
    assert paths["fasta"].is_relative_to(store / snapshot["snapshot_path"] / "source")
    assert set(_hashes(prepared)) == {
        "sequences.fai",
        "sequences.tsv",
        "provenance.json",
        "README.md",
    }, "The prepared view must not duplicate or rewrite the native FASTA"
    for field in ("provenance", "readme"):
        assert str(store) not in paths[field].read_text()
    return prepared


def _move_and_read(
    tmp_path: Path,
    source: Path,
    store: Path,
    snapshot: dict[str, Any],
    prepared: Path,
) -> dict[str, Any]:
    """Return an independent summary with no old source/store, database or BioV."""
    python = _configured_path("BIOV_STANDALONE_FASTA", executable=True)
    snapshot_relative = snapshot["snapshot_path"]
    prepared_relative = prepared.relative_to(store)
    expected_source = _hashes(source)
    expected_snapshot = _hashes(store / snapshot_relative)
    expected_prepared = _hashes(prepared)
    with tempfile.TemporaryDirectory(prefix="independent-prepared-fasta-") as directory:
        consumer = Path(directory).resolve()
        assert not consumer.is_relative_to(_ROOT)
        assert not consumer.is_relative_to(tmp_path)
        moved = consumer / "portable-store"
        shutil.copytree(store / snapshot_relative, moved / snapshot_relative)
        shutil.copytree(prepared, moved / prepared_relative)
        assert not (moved / "catalog").exists()
        assert not list(moved.rglob("*.sqlite*"))
        reader = consumer / "read_prepared_fasta.py"
        shutil.copyfile(_ROOT / "tests/_prepared_fasta_reader.py", reader)
        # These are disposable test copies, never the user's original package.
        assert _hashes(source) == expected_source
        shutil.rmtree(source)
        shutil.rmtree(store)
        assert not source.exists() and not store.exists()
        before = _hashes(moved)
        result = _run_offline(
            [str(python), "-I", str(reader), str(moved / prepared_relative)], consumer
        )
        assert _hashes(moved) == before, (
            "Independent access must not rebuild/write indexes"
        )
        assert _hashes(moved / snapshot_relative) == expected_snapshot
        assert _hashes(moved / prepared_relative) == expected_prepared
        summary = json.loads(result.stdout)
        assert summary["network_blocked"] is True
        assert summary["inventory_verified"] is True
        assert summary["readme_example_executed"] is True
        return summary


@pytest.mark.parametrize("newline", [b"\n", b"\r\n"], ids=["lf", "crlf"])
def test_prepared_fasta_synthetic_offline_move(tmp_path: Path, newline: bytes) -> None:
    """Reuse and move case-preserving wrapped indexes into an independent reader."""
    binary = _binary()
    source = tmp_path / "imports/native-package"
    _synthetic_source(source, newline)
    source_before = _hashes(source)
    store = tmp_path / "producer-store"
    store.mkdir()
    snapshot = _register(binary, tmp_path, store, source)
    saved_snapshot = _hashes(store / snapshot["snapshot_path"])
    request = _request(snapshot)
    result = _cli(binary, tmp_path, store, request)
    assert result is not None and result["reused"] is False
    prepared = _check_result(result, store, snapshot)
    assert result["sequence_count"] == len(_SEQUENCES)
    assert result["total_bases"] == sum(map(len, _SEQUENCES.values()))
    saved_prepared = _hashes(prepared)
    repeated = _cli(binary, tmp_path, store, request)
    assert repeated == dict(result, reused=True)
    assert _hashes(prepared) == saved_prepared
    assert _hashes(store / snapshot["snapshot_path"]) == saved_snapshot
    assert _hashes(source) == source_before
    summary = _move_and_read(tmp_path, source, store, snapshot, prepared)
    assert summary["sequences"] == [
        {
            "id": name,
            "length": len(sequence),
            "sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        }
        for name, sequence in _SEQUENCES.items()
    ]
    assert summary["total_bases"] == sum(map(len, _SEQUENCES.values()))
    assert summary["masked_lowercase_bases"] == sum(
        sum(base.islower() for base in sequence) for sequence in _SEQUENCES.values()
    )
    assert summary["slice_checks"] >= len(_SEQUENCES) * 4


@pytest.mark.parametrize(
    "target", ["fai", "dictionary", "provenance", "readme", "fasta"]
)
def test_prepared_fasta_rejects_altered_outputs_or_source(
    tmp_path: Path, target: str
) -> None:
    """Never reuse or silently replace corrupted output or immutable input bytes."""
    binary = _binary()
    source = tmp_path / "imports/native-package"
    _synthetic_source(source)
    originals = _hashes(source)
    store = tmp_path / "producer-store"
    store.mkdir()
    snapshot = _register(binary, tmp_path, store, source)
    request = _request(snapshot)
    result = _cli(binary, tmp_path, store, request)
    assert result is not None
    _check_result(result, store, snapshot)
    path = _relative_path(store, result[target]["relative_path"])
    original = path.read_bytes()
    if target == "provenance":
        record = json.loads(original)
        record["sequence_count"] += 1
        path.write_text(json.dumps(record), encoding="utf-8")
    elif target == "fasta":
        path.write_bytes(original.replace(b"ACGTacg", b"TCGTacg", 1))
    else:
        path.write_bytes(original + b"altered output\n")
    assert path.read_bytes() != original
    corrupted = _hashes(store)
    _cli(binary, tmp_path, store, request, success=False)
    assert _hashes(store) == corrupted, (
        "Rejection must preserve the evidence/prior snapshot"
    )
    assert _hashes(source) == originals


def test_prepared_fasta_changed_snapshot_changes_recipe(tmp_path: Path) -> None:
    """Exact snapshot identity participates even when selected FASTA bytes match."""
    binary = _binary()
    source = tmp_path / "imports/native-package"
    _synthetic_source(source)
    store = tmp_path / "producer-store"
    store.mkdir()
    first_snapshot = _register(binary, tmp_path, store, source)
    first = _cli(binary, tmp_path, store, _request(first_snapshot))
    assert first is not None
    first_prepared = _check_result(first, store, first_snapshot)
    first_before = _hashes(first_prepared)
    snapshot_before = _hashes(store / first_snapshot["snapshot_path"])
    with (source / "README.md").open("a", encoding="utf-8") as stream:
        stream.write("New native snapshot; selected FASTA remains identical.\n")
    second_snapshot = _register(binary, tmp_path, store, source)
    second = _cli(binary, tmp_path, store, _request(second_snapshot))
    assert second is not None and second["reused"] is False
    _check_result(second, store, second_snapshot)
    assert first_snapshot["snapshot_id"] != second_snapshot["snapshot_id"]
    assert first["recipe_id"] != second["recipe_id"]
    assert _hashes(first_prepared) == first_before
    assert _hashes(store / first_snapshot["snapshot_path"]) == snapshot_before
    assert _cli(binary, tmp_path, store, _request(first_snapshot)) == dict(
        first, reused=True
    )


def test_prepared_fasta_real_refseq_offline_move(tmp_path: Path) -> None:
    """Read the whole actual RefSeq genome and slices without changing originals."""
    if "BIOV_TEST_REFSEQ_PACKAGE" not in os.environ:
        pytest.skip("Set BIOV_TEST_REFSEQ_PACKAGE for the existing real RefSeq package")
    binary = _binary()
    original = _configured_path("BIOV_TEST_REFSEQ_PACKAGE")
    before = _hashes(original)
    source = tmp_path / "imports/real-refseq"
    shutil.copytree(original, source)
    try:
        store = tmp_path / "producer-store"
        store.mkdir()
        snapshot = _register(binary, tmp_path, store, source)
        source_path = (
            f"ncbi_dataset/data/{_ACCESSION}/{_ACCESSION}_ASM584v2_genomic.fna"
        )
        result = _cli(binary, tmp_path, store, _request(snapshot, source_path))
        assert result is not None and result["reused"] is False
        prepared = _check_result(result, store, snapshot)
        assert result["sequence_count"] == 1 and result["total_bases"] == 4641652
        summary = _move_and_read(tmp_path, source, store, snapshot, prepared)
        assert [
            (sequence["id"], sequence["length"]) for sequence in summary["sequences"]
        ] == [("NC_000913.3", 4641652)]
        assert summary["total_bases"] == 4641652
    finally:
        assert _hashes(original) == before, "Original downloaded RefSeq package changed"


def test_prepared_fasta_interrupted_stage_retries(tmp_path: Path) -> None:
    """SIGKILL during indexing leaves an ignored stage and a retryable snapshot."""
    binary = _binary()
    source = tmp_path / "imports/interruption-package"
    fasta = source / _FASTA_PATH
    fasta.parent.mkdir(parents=True)
    (source / "README.md").write_text(
        "Generated interruption fixture, not provider data.\n"
    )
    line = b"ACGTacgtnN" * 8 + b"\n"
    line_count = (64 * 1024 * 1024) // len(line) + 1
    with fasta.open("wb") as stream:
        stream.write(b">interruption.1\n")
        block = line * 1024
        complete, remainder = divmod(line_count, 1024)
        for _ in range(complete):
            stream.write(block)
        stream.write(line * remainder)
    _native_metadata(source, fasta)
    originals = _hashes(source)
    store = tmp_path / "producer-store"
    store.mkdir()
    snapshot = _register(binary, tmp_path, store, source)
    snapshot_before = _hashes(store / snapshot["snapshot_path"])
    request = _request(snapshot)
    request_path = tmp_path / "preparation.json"
    request_path.write_text(json.dumps(request), encoding="utf-8")
    process = subprocess.Popen(
        _offline_command(
            [
                str(binary),
                "prepared",
                "fasta",
                "--store-root",
                str(store),
                "--request-file",
                str(request_path),
            ]
        ),
        cwd=tmp_path,
        env=_environment(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    interrupted = False
    deadline = time.monotonic() + 30
    try:
        while time.monotonic() < deadline and process.poll() is None:
            if any((store / ".staging").glob("prepare-fasta-*")):
                process.kill()
                interrupted = True
                break
            time.sleep(0.002)
        stdout, stderr = process.communicate(timeout=10)
        assert interrupted and process.returncode < 0, (
            "Could not interrupt inside the prepared staging boundary",
            stdout,
            stderr,
        )
    finally:
        if process.poll() is None:
            process.kill()
            process.communicate(timeout=10)
    stale = sorted((store / ".staging").glob("prepare-fasta-*"))
    assert stale, "A real SIGKILL should leave a stage for retry acceptance"
    stale_before = {path.name: _hashes(path) for path in stale}
    assert not list((store / "prepared").glob("sha256-*")), "No incomplete READY view"
    assert _hashes(store / snapshot["snapshot_path"]) == snapshot_before
    completed = _cli(binary, tmp_path, store, request)
    assert completed is not None and completed["reused"] is False
    _check_result(completed, store, snapshot)
    assert completed["sequence_count"] == 1
    assert completed["total_bases"] == line_count * 80
    assert {path.name: _hashes(path) for path in stale} == stale_before
    assert _hashes(store / snapshot["snapshot_path"]) == snapshot_before
    assert _hashes(source) == originals
    assert _cli(binary, tmp_path, store, request) == dict(completed, reused=True)
