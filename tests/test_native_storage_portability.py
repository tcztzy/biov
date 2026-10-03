"""Installed CLI acceptance for portable native storage, with no network downloads.

CI enables the synthetic test with BIOV_TEST_BINARY and BIOV_STANDALONE_PYTHON
(a separate pyarrow-only virtualenv, shared with Arrow portability acceptance).
Real package tests are opt-in: BIOV_TEST_REFSEQ_PACKAGE names the extracted
GCF_000005845.2 package; BIOV_TEST_PDB_PACKAGE names a directory containing
1crn.cif.gz, 1crn.cif and 1crn.entry.json. The latter additionally requires
BIOV_STANDALONE_BIOPYTHON (only biopython==1.88 and its numpy dependency).

Original fixture trees are never modified. Only disposable copies are registered.
Configured-but-broken paths/dependencies are failures, not skipped acceptance.
The kernel network sandbox needs Linux and libseccomp.so.2 (Ubuntu CI provides
it); no host settings, privileges or network namespaces are changed. Both the
Rust producer/resolver and standalone Python consumers inherit the restriction.
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
_TIMEOUT = 90
_REFSEQ = "refseq.gcf:GCF_000005845.2"
_PDB_FILES = ("1crn.cif", "1crn.cif.gz", "1crn.entry.json")
_BASE_ENV = ("BIOV_TEST_BINARY", "BIOV_STANDALONE_PYTHON")


def _configured_path(name: str, *, executable: bool = False) -> Path:
    """Return a valid explicit path, preserving virtualenv executable symlinks."""
    value = os.environ.get(name)
    assert value, f"Required installed acceptance configuration is missing: {name}"
    path = Path(value).absolute()
    if executable:
        assert path.is_file() and os.access(path, os.X_OK), (
            f"{name} is not executable: {path}"
        )
    else:
        assert path.is_dir(), f"{name} is not an existing fixture directory: {path}"
    return path


def _run_offline(arguments: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    """Return one bounded child result under kernel-enforced network denial."""
    assert sys.platform == "linux", (
        "Offline acceptance requires Linux + libseccomp.so.2"
    )
    try:
        ctypes.CDLL("libseccomp.so.2")
    except OSError as error:
        pytest.fail(f"Offline acceptance needs libseccomp.so.2: {error}")
    result = subprocess.run(
        [
            sys.executable,
            "-I",
            str(_ROOT / "tests/_native_storage_sandbox.py"),
            *arguments,
        ],
        cwd=cwd,
        env={
            "PATH": "",
            "PYTHONPATH": "",
            "PYTHONNOUSERSITE": "1",
            "POLARS_MAX_THREADS": "2",
        },
        capture_output=True,
        text=True,
        timeout=_TIMEOUT,
        check=False,
    )
    assert result.returncode == 0, (
        f"Offline acceptance command failed: {arguments[:3]}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    return result


def _cli(
    binary: Path,
    cwd: Path,
    operation: str,
    store: Path,
    request: dict[str, Any],
    source_root: Path | None = None,
) -> dict[str, Any]:
    """Return CLI JSON without importing BioV or duplicating an MCP client."""
    with tempfile.TemporaryDirectory(prefix="storage-request-", dir=cwd) as directory:
        request_path = Path(directory) / "request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        arguments = [str(binary), "storage", operation, "--store-root", str(store)]
        if source_root is not None:
            arguments += ["--source-root", str(source_root)]
        arguments += ["--request-file", str(request_path)]
        return json.loads(_run_offline(arguments, cwd).stdout)


def _file_hashes(root: Path) -> dict[str, tuple[int, str]]:
    """Return streamed identities for every regular file in a tree."""
    output = {}
    for path in sorted(root.rglob("*")):
        assert not path.is_symlink(), (
            f"Acceptance fixture must be self-contained: {path}"
        )
        if path.is_file():
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            output[path.relative_to(root).as_posix()] = (path.stat().st_size, digest)
    return output


def _portable_path(root: Path, name: str) -> Path:
    """Return an existing confined file or directory from a store-relative path."""
    assert name and not Path(name).is_absolute() and "\\" not in name, name
    assert all(part not in {"", ".", ".."} for part in name.split("/")), name
    path = root / name
    assert path.exists() and path.resolve().is_relative_to(root.resolve()), name
    return path


def _ready_paths(result: dict[str, Any], store: Path) -> list[Path]:
    """Check that resolved outputs are directly readable ordinary local paths.

    Returns:
        The complete selected representation's file paths on the current host.
    """
    assert result["status"] == "ready", result
    paths = []
    for item in result["paths"]:
        path = _portable_path(store, item["relative_path"])
        assert path.is_file() and not path.is_symlink()
        assert Path(item["execution_host_path"]) == path
        with path.open("rb") as stream:
            assert stream.read(1), "Resolved file must contain real data"
        paths.append(path)
    assert paths
    return paths


def _synthetic_refseq(root: Path) -> None:
    """Write a tiny materialized native-layout fixture without provider downloads."""
    data = root / "ncbi_dataset/data"
    accession = data / "GCF_000005845.2"
    accession.mkdir(parents=True)
    (root / "empty-native-directory").mkdir()
    (root / "README.md").write_text(
        "Synthetic NCBI-layout fixture, not provider data.\n"
    )
    (root / "μ-native-note.txt").write_text("Unicode path inventory fixture.\n")
    fasta = accession / "genome.fna"
    fasta.write_text(
        ">synthetic_a example\nACGTNN\n>synthetic_b\nGGCC\n", encoding="utf-8"
    )
    catalog = data / "dataset_catalog.json"
    catalog.write_text(
        json.dumps(
            {
                "apiVersion": "V2",
                "assemblies": [
                    {
                        "accession": "GCF_000005845.2",
                        "files": [
                            {
                                "filePath": "GCF_000005845.2/genome.fna",
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
    (root / "md5sum.txt").write_text(
        "".join(
            hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest()
            + "  "
            + path.relative_to(root).as_posix()
            + "\n"
            for path in (catalog, fasta)
        ),
        encoding="utf-8",
    )


def _exercise_portability(
    tmp_path: Path,
    source: Path,
    mode: str,
    request: dict[str, Any],
    representation: str,
    *,
    create_revision: bool = False,
) -> dict[str, Any]:
    """Register, relocate, rediscover and independently read an immutable copy.

    Returns:
        The standalone consumer's verified full-data scientific summary.
    """
    binary = _configured_path("BIOV_TEST_BINARY", executable=True)
    python = _configured_path(
        "BIOV_STANDALONE_BIOPYTHON" if mode == "pdb" else "BIOV_STANDALONE_PYTHON",
        executable=True,
    )
    original_hashes = _file_hashes(source)
    original_store = tmp_path / "producer-store"
    original_store.mkdir()
    registration = _cli(
        binary, tmp_path, "register", original_store, request, source.parent
    )
    assert registration["reused"] is False
    snapshot = registration["snapshot"]
    assert snapshot["file_count"] == len(original_hashes)
    assert snapshot["total_bytes"] == sum(size for size, _ in original_hashes.values())
    original_snapshot = _portable_path(original_store, snapshot["snapshot_path"])
    receipt_path = _portable_path(original_store, snapshot["receipt_path"])
    assert receipt_path == original_snapshot / "acquisition.json"
    assert _file_hashes(original_snapshot / "source") == original_hashes
    saved_envelope = _file_hashes(original_snapshot)
    for companion in ("README.md", "acquisition.json", "checksums.sha256"):
        text = (original_snapshot / companion).read_text()
        assert str(source) not in text and str(original_store) not in text
    assert (
        _cli(binary, tmp_path, "register", original_store, request, source.parent)[
            "reused"
        ]
        is True
    )
    assert _file_hashes(original_snapshot) == saved_envelope
    resolve = {
        "reference": request["canonical_ref"],
        "representation": representation,
        "snapshot_id": None,
        "scope": None,
    }
    first = _cli(binary, tmp_path, "resolve", original_store, resolve)
    _ready_paths(first, original_store)
    assert first["snapshot"] == snapshot
    miss = dict(resolve, reference="pdb:9ZZZ")
    assert _cli(binary, tmp_path, "resolve", original_store, miss) == {"status": "miss"}
    if mode in {"synthetic", "refseq"}:
        missing_rna = dict(resolve, representation="rna_fasta")
        unavailable = _cli(binary, tmp_path, "resolve", original_store, missing_rna)
        assert unavailable["status"] == "unavailable"
        assert "rna_fasta" not in unavailable["available_representations"]
    if create_revision:
        with (source / "README.md").open("a", encoding="utf-8") as stream:
            stream.write("Second immutable synthetic native tree.\n")
        changed = _cli(
            binary, tmp_path, "register", original_store, request, source.parent
        )
        assert changed["snapshot"]["snapshot_id"] != snapshot["snapshot_id"]
        assert _file_hashes(original_snapshot) == saved_envelope
        ambiguous = _cli(binary, tmp_path, "resolve", original_store, resolve)
        assert ambiguous["status"] == "ambiguous"
        assert {item["snapshot_id"] for item in ambiguous["candidates"]} == {
            snapshot["snapshot_id"],
            changed["snapshot"]["snapshot_id"],
        }
    resolve["snapshot_id"] = snapshot["snapshot_id"]

    # This root is unrelated to both the checkout and producer temp directory.
    # No BioV database, response cache or historical source location is needed.
    with tempfile.TemporaryDirectory(
        prefix="independent-native-consumer-"
    ) as directory:
        consumer = Path(directory).resolve()
        assert not consumer.is_relative_to(_ROOT)
        assert not consumer.is_relative_to(tmp_path)
        moved_store = consumer / "moved-store"
        shutil.move(original_store, moved_store)
        shutil.rmtree(source)
        assert not source.exists() and not original_store.exists()
        catalog = moved_store / "catalog"
        if catalog.exists():
            shutil.rmtree(catalog)
        assert not list(moved_store.rglob("*.sqlite*")), "Reader receives no database"
        # A new CLI process reconstructs discovery from portable receipts only.
        moved = _cli(binary, consumer, "resolve", moved_store, resolve)
        paths = _ready_paths(moved, moved_store)
        assert moved["snapshot"] == snapshot
        assert all(not str(path).startswith(str(original_store)) for path in paths)
        moved_snapshot = _portable_path(moved_store, snapshot["snapshot_path"])
        assert _file_hashes(moved_snapshot) == saved_envelope
        if create_revision:
            ambiguous = _cli(
                binary,
                consumer,
                "resolve",
                moved_store,
                dict(resolve, snapshot_id=None),
            )
            assert ambiguous["status"] == "ambiguous"
        # Also copy just one self-contained envelope; its ordinary reader gets
        # no store, registry, CLI response, source root, or original file names.
        standalone = consumer / "exported-snapshot" / moved_snapshot.name
        shutil.copytree(moved_snapshot, standalone)
        reader = consumer / "read_native_snapshot.py"
        shutil.copyfile(_ROOT / "tests/_native_storage_reader.py", reader)
        arguments = [str(python), "-I", str(reader), mode, str(standalone)]
        if mode == "refseq":
            refseq_reader = consumer / "inspect_refseq_example.py"
            shutil.copyfile(_ROOT / "scripts/inspect_refseq_example.py", refseq_reader)
            arguments += ["--refseq-reader", str(refseq_reader)]
        shutil.rmtree(moved_store)
        assert not moved_store.exists(), (
            "Standalone consumer must not depend on the store"
        )
        result = json.loads(_run_offline(arguments, consumer).stdout)
        assert result["inventory_verified"] is True
        assert result["network_blocked"] is True
        assert _file_hashes(standalone) == saved_envelope
        return result["summary"]


def test_native_storage_synthetic_offline_move(tmp_path: Path) -> None:
    """CI verifies offline relocation, immutable revisions and explicit selection."""
    if not any(name in os.environ for name in _BASE_ENV):
        pytest.skip(
            "Set BIOV_TEST_BINARY and BIOV_STANDALONE_PYTHON for installed acceptance"
        )
    source = tmp_path / "imports" / "synthetic-native"
    _synthetic_refseq(source)
    summary = _exercise_portability(
        tmp_path,
        source,
        "synthetic",
        {
            "source_path": source.name,
            "requested_ref": _REFSEQ,
            "canonical_ref": _REFSEQ,
            "declaration": {"provider": "refseq"},
        },
        "genome_fasta",
        create_revision=True,
    )
    assert summary == {"fasta_records": 2, "total_bases": 10, "unambiguous_records": 1}


def test_native_storage_real_refseq_offline_move(tmp_path: Path) -> None:
    """Preserve and analyze every actual RefSeq file without modifying the fixture."""
    if "BIOV_TEST_REFSEQ_PACKAGE" not in os.environ:
        pytest.skip(
            "Set BIOV_TEST_REFSEQ_PACKAGE to opt into the existing real package"
        )
    original = _configured_path("BIOV_TEST_REFSEQ_PACKAGE")
    before = _file_hashes(original)
    assert len(before) == 9, (
        "Expected the full inspected GCF_000005845.2 native package"
    )
    source = tmp_path / "imports" / "real-refseq-native"
    try:
        shutil.copytree(original, source)
        summary = _exercise_portability(
            tmp_path,
            source,
            "refseq",
            {
                "source_path": source.name,
                "requested_ref": _REFSEQ,
                "canonical_ref": _REFSEQ,
                "declaration": {"provider": "refseq"},
            },
            "genome_fasta",
        )
        assert summary["fasta_records"] == {"genome": 1, "protein": 4300, "cds": 4318}
    finally:
        assert _file_hashes(original) == before, "Original real package was modified"


def test_native_storage_real_pdb_offline_move(tmp_path: Path) -> None:
    """Use stock Biopython on real native PDB files, excluding investigation tools."""
    if "BIOV_TEST_PDB_PACKAGE" not in os.environ:
        pytest.skip("Set BIOV_TEST_PDB_PACKAGE to opt into the existing real PDB files")
    original = _configured_path("BIOV_TEST_PDB_PACKAGE")
    before = _file_hashes(original)
    assert set(_PDB_FILES).issubset(before), f"PDB fixture must contain {_PDB_FILES}"
    source = tmp_path / "imports" / "real-pdb-native"
    source.mkdir(parents=True)
    try:
        for name in _PDB_FILES:
            shutil.copyfile(original / name, source / name)
        summary = _exercise_portability(
            tmp_path,
            source,
            "pdb",
            {
                "source_path": source.name,
                "requested_ref": "pdb:1CRN",
                "canonical_ref": "pdb:1CRN",
                "declaration": {
                    "provider": "pdb",
                    "scope": "entry",
                    "representations": {
                        "structure_cif": ["1crn.cif"],
                        "structure_cif_gzip": ["1crn.cif.gz"],
                        "entry_json": ["1crn.entry.json"],
                    },
                },
            },
            "structure_cif",
        )
        assert summary["models"] == 1 and summary["residues"] == 46
        assert summary["atoms"] == 327
    finally:
        assert _file_hashes(original) == before, (
            "Original PDB inspection tree was modified"
        )


def test_native_storage_interrupted_registration_retries(tmp_path: Path) -> None:
    """A killed staged copy stays a miss; retry handles a native file over 64 MiB."""
    if not any(name in os.environ for name in _BASE_ENV):
        pytest.skip("Set BIOV_TEST_BINARY for installed native storage acceptance")
    binary = _configured_path("BIOV_TEST_BINARY", executable=True)
    assert sys.platform == "linux", "Crash acceptance requires Linux + libseccomp.so.2"
    ctypes.CDLL("libseccomp.so.2")
    source_root = tmp_path / "imports"
    source = source_root / "synthetic-native"
    source.mkdir(parents=True)
    payload = source / "synthetic.bin"
    # Deliberately opaque synthetic bytes. The PDB route verifies caller-declared
    # paths and local integrity only; this is no claim that these bytes are PDB.
    size = 64 * 1024 * 1024 + 1
    with payload.open("wb") as stream:
        stream.write(b"synthetic interrupted copy fixture\n")
        stream.truncate(size)
    before = _file_hashes(source)
    store = tmp_path / "native-store"
    store.mkdir()
    request = {
        "source_path": source.name,
        "requested_ref": "pdb:1ABC",
        "canonical_ref": "pdb:1ABC",
        "declaration": {
            "provider": "pdb",
            "scope": "entry",
            "representations": {"opaque_synthetic_fixture": [payload.name]},
        },
    }
    request_path = tmp_path / "registration.json"
    request_path.write_text(json.dumps(request), encoding="utf-8")
    process = subprocess.Popen(
        [
            sys.executable,
            "-I",
            str(_ROOT / "tests/_native_storage_sandbox.py"),
            str(binary),
            "storage",
            "register",
            "--store-root",
            str(store),
            "--source-root",
            str(source_root),
            "--request-file",
            str(request_path),
        ],
        cwd=tmp_path,
        env={"PATH": "", "PYTHONPATH": "", "POLARS_MAX_THREADS": "2"},
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    interrupted = False
    deadline = time.monotonic() + 20
    try:
        while time.monotonic() < deadline and process.poll() is None:
            if any((store / ".staging").glob("register-*/source/synthetic.bin")):
                process.kill()
                interrupted = True
                break
            time.sleep(0.002)
        stdout, stderr = process.communicate(timeout=10)
        assert interrupted and process.returncode < 0, (
            "Could not interrupt registration inside its staging boundary",
            stdout,
            stderr,
        )
    finally:
        if process.poll() is None:
            process.kill()
            process.communicate(timeout=10)
    staging = store / ".staging"
    leftovers = _file_hashes(staging)
    assert leftovers, "SIGKILL must leave incomplete staging for recovery acceptance"
    assert _file_hashes(source) == before
    resolution = {
        "reference": "pdb:1ABC",
        "representation": "opaque_synthetic_fixture",
        "snapshot_id": None,
        "scope": "entry",
    }
    assert _cli(binary, tmp_path, "resolve", store, resolution) == {"status": "miss"}
    assert not list((store / "artifacts").rglob("acquisition.json"))
    completed = _cli(binary, tmp_path, "register", store, request, source_root)
    assert completed["reused"] is False
    assert completed["snapshot"]["total_bytes"] == size > 64 * 1024 * 1024
    paths = _ready_paths(_cli(binary, tmp_path, "resolve", store, resolution), store)
    assert len(paths) == 1 and paths[0].stat().st_size == size
    assert _file_hashes(paths[0].parent) == before
    assert _file_hashes(source) == before
    assert _file_hashes(staging) == leftovers, (
        "Recovery must not treat/delete stale work as READY"
    )
