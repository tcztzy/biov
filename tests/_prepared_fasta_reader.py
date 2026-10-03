"""Independent prepared-FASTA acceptance reader copied outside the BioV checkout.

Biopython performs the sequential reference parse. Unmodified pysam/HTSlib
performs indexed access through the explicit, externally located FAI. This
helper verifies one documented portable export; it is not a BioV runtime API or
a replacement FASTA parser. All reads are offline and must leave files intact.
"""

import argparse
import ast
import contextlib
import csv
import errno
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import os
import random
import re
import socket
import struct
import sys
from pathlib import Path
from typing import Any

_EXPECTED_PACKAGES = {"biopython": "1.88", "numpy": "2.5.3", "pysam": "0.24.1"}


def _reader_environment() -> None:
    """Verify isolated packages and OS-level denial of native network access.

    Raises:
        AssertionError: If isolation or the inherited kernel filter is missing.
    """
    assert os.environ["PATH"] == os.environ["PYTHONPATH"] == ""
    assert sys.prefix != sys.base_prefix, "Reader needs its own clean virtualenv"
    assert importlib.util.find_spec("biov") is None, "BioV must not be installed"
    installed = {
        distribution.metadata["Name"].lower(): distribution.version
        for distribution in importlib.metadata.distributions()
    }
    assert installed == _EXPECTED_PACKAGES, (installed, _EXPECTED_PACKAGES)
    for family, kind in (
        (socket.AF_INET, socket.SOCK_STREAM),
        (socket.AF_INET6, socket.SOCK_DGRAM),
        (socket.AF_UNIX, socket.SOCK_STREAM),
    ):
        try:
            connection = socket.socket(family, kind)
        except PermissionError as error:
            assert error.errno == errno.EPERM
        else:
            connection.close()
            raise AssertionError("The inherited kernel network filter is missing")

    def independent_only(event: str, args: tuple[Any, ...]) -> None:
        """Reject hidden service/database dependencies and escaped child execution.

        Raises:
            AssertionError: If a forbidden operation or BioV import is attempted.
        """
        if event.startswith(("socket.", "subprocess.", "sqlite3.")) or event in {
            "os.system",
            "os.exec",
            "os.posix_spawn",
            "os.fork",
            "os.forkpty",
        }:
            raise AssertionError(f"Reader must remain independent and offline: {event}")
        if event == "import" and args[0].split(".")[0] == "biov":
            raise AssertionError(f"Forbidden reader dependency: {args[0]}")

    sys.addaudithook(independent_only)


def _digest(path: Path) -> str:
    """Return SHA-256 without retaining a file's payload."""
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _safe_file(root: Path, name: str) -> Path:
    """Return a non-symlink ordinary file with a portable relative name."""
    assert name and not Path(name).is_absolute() and "\\" not in name, name
    assert all(part not in {"", ".", ".."} for part in name.split("/")), name
    path = root / name
    assert path.is_file() and not path.is_symlink(), name
    assert path.resolve().is_relative_to(root.resolve()), name
    return path


def _native_inventory(snapshot: Path) -> dict[str, Any]:
    """Verify every native file and directory, including the canonical tree hash.

    Returns:
        The ordinary source receipt after independent byte-identity verification.
    """
    receipt = json.loads(_safe_file(snapshot, "acquisition.json").read_text())
    assert receipt["schema_version"] == 1 and receipt["origin"] == "local_copy"
    source = snapshot / "source"
    inventory = receipt["inventory"]
    names = [entry["path"] for entry in inventory]
    assert names == sorted(set(names), key=lambda name: name.encode("utf-8"))
    assert set(names) == {
        path.relative_to(source).as_posix() for path in source.rglob("*")
    }
    tree_hash = hashlib.sha256(b"biov-native-tree-v1\0")
    checksums = {}
    for entry in inventory:
        name = entry["path"]
        assert not (source / name).is_symlink(), name
        raw_name = name.encode("utf-8")
        if entry["kind"] == "directory":
            assert (source / name).is_dir()
            assert entry["bytes"] == 0 and entry["sha256"] is None
            kind, size, digest = b"D", 0, bytes(32)
        else:
            assert entry["kind"] == "file", entry
            path = _safe_file(source, name)
            sha256 = _digest(path)
            size = path.stat().st_size
            assert entry["bytes"] == size and entry["sha256"] == sha256, name
            kind, digest = b"F", bytes.fromhex(sha256)
            checksums["source/" + name] = sha256
        tree_hash.update(kind + struct.pack(">I", len(raw_name)) + raw_name)
        tree_hash.update(struct.pack(">Q", size) + digest)
    assert receipt["source_content_sha256"] == tree_hash.hexdigest()
    assert receipt["snapshot_id"] == "sha256-" + tree_hash.hexdigest() == snapshot.name
    observed = {}
    for line in _safe_file(snapshot, "checksums.sha256").read_text().splitlines():
        checksum_digest, name = line.split(maxsplit=1)
        assert name not in observed
        observed[name] = checksum_digest
    assert observed == checksums
    assert _digest(_safe_file(snapshot, "README.md")) == receipt["readme_sha256"]
    return receipt


def _prepared_inventory(prepared: Path) -> tuple[dict[str, Any], Path]:
    """Validate the documented relative dependency closure and output identities.

    Returns:
        The prepared record and exact FASTA selected from its immutable snapshot.
    """
    store = prepared.parent.parent
    assert prepared.parent.name == "prepared"
    assert not list(store.rglob("*.sqlite*")) and not (store / "catalog").exists()
    record = json.loads(_safe_file(prepared, "provenance.json").read_text())
    assert record["schema_version"] == 1
    assert record["recipe_id"] == prepared.name
    recipe = record["recipe"]
    canonical_recipe = json.dumps(
        recipe, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    recipe_digest = hashlib.sha256(
        b"biov-prepared-fasta-recipe-v1\0" + canonical_recipe
    ).hexdigest()
    assert record["recipe_id"] == "sha256-" + recipe_digest
    assert recipe["schema_version"] == 1
    assert recipe["reference"] == "refseq.gcf:GCF_000005845.2"
    assert recipe["snapshot_id"].startswith("sha256-")
    assert recipe["source_bytes"] > 0
    assert recipe["implementation"] and recipe["implementation_version"]
    assert recipe["library"] == "noodles-fasta" and recipe["library_version"]
    # This protocol permits ../../ solely for these explicit closure dependencies.
    for key in ("source_relative_path", "source_snapshot_relative_path"):
        relative = record[key]
        assert relative.startswith("../../artifacts/") and "\\" not in relative
        assert not Path(relative).is_absolute()
    snapshot = (prepared / record["source_snapshot_relative_path"]).resolve()
    fasta = (prepared / record["source_relative_path"]).resolve()
    assert snapshot.is_relative_to(store / "artifacts") and snapshot.is_dir()
    assert fasta == _safe_file(snapshot / "source", recipe["source_path"]).resolve()
    receipt = _native_inventory(snapshot)
    assert receipt["snapshot_id"] == recipe["snapshot_id"]
    assert receipt["canonical_ref"] == recipe["reference"]
    assert receipt["source_content_sha256"] == recipe["source_content_sha256"]
    assert recipe["source_path"] in receipt["representations"]["genome_fasta"]
    assert _digest(fasta) == recipe["source_sha256"]
    assert fasta.stat().st_size == recipe["source_bytes"]
    outputs = record["outputs"]
    assert {entry["path"] for entry in outputs} == {
        "sequences.fai",
        "sequences.tsv",
        "README.md",
    }
    assert len(outputs) == 3
    for entry in outputs:
        path = _safe_file(prepared, entry["path"])
        assert path.stat().st_size == entry["bytes"]
        assert _digest(path) == entry["sha256"]
    assert {path.name for path in prepared.iterdir()} == {
        "sequences.fai",
        "sequences.tsv",
        "README.md",
        "provenance.json",
    }
    assert not Path(str(fasta) + ".fai").exists(), (
        "Index must remain an external derived view"
    )
    return record, fasta


def _run_readme(prepared: Path) -> dict[str, Any]:
    """Return the summary from the emitted ordinary-reader example."""
    examples = re.findall(
        r"^```python\s*\n(.*?)^```\s*$",
        (prepared / "README.md").read_text(),
        re.MULTILINE | re.DOTALL,
    )
    assert len(examples) == 1, (
        "Prepared README needs one executable ordinary-reader example"
    )
    previous = Path.cwd()
    output = io.StringIO()
    try:
        os.chdir(prepared)
        with contextlib.redirect_stdout(output):
            exec(  # noqa: S102 - Trusted BioV-emitted example, never provider source code.
                compile(examples[0], "prepared README example", "exec"),
                {"__name__": "__main__"},
            )
    finally:
        os.chdir(previous)
    result = ast.literal_eval(output.getvalue().strip())
    assert isinstance(result, dict), "Example must report a useful sequence summary"
    return result


def _read_sequences(
    prepared: Path, fasta: Path, record: dict[str, Any]
) -> dict[str, Any]:
    """Compare complete independent sequential records against indexed access.

    Returns:
        Record identities, boundary/interior fetch coverage and a full-data summary.
    """
    import pysam
    from Bio import SeqIO

    # Established independent scientific parsers, not a hand-written FASTA decoder.
    sequential = list(SeqIO.parse(fasta, "fasta"))
    ids = [sequence.id for sequence in sequential]
    lengths = [len(sequence) for sequence in sequential]
    assert len(set(ids)) == len(ids) and sequential
    with _safe_file(prepared, "sequences.tsv").open(newline="") as stream:
        table = csv.DictReader(stream, delimiter="\t")
        assert table.fieldnames == ["sequence_id", "length"]
        dictionary = list(table)
    assert dictionary == [
        {"sequence_id": name, "length": str(length)}
        for name, length in zip(ids, lengths, strict=True)
    ]
    with _safe_file(prepared, "sequences.fai").open(newline="") as stream:
        fai = list(csv.reader(stream, delimiter="\t"))
    assert len(fai) == len(ids)
    for row, name, length in zip(fai, ids, lengths, strict=True):
        assert len(row) == 5 and row[0] == name and int(row[1]) == length
        assert int(row[2]) > 0 and int(row[3]) > 0 and int(row[4]) >= int(row[3])
    assert record["sequence_count"] == len(sequential)
    assert record["total_bases"] == sum(lengths)
    summaries = []
    slice_checks = total_gc = lowercase = ambiguous = 0
    generator = random.Random(0)  # noqa: S311 - Reproducible test coordinates, not security.
    with pysam.FastaFile(
        str(fasta), filepath_index=str(prepared / "sequences.fai")
    ) as indexed:
        assert list(indexed.references) == ids
        assert list(indexed.lengths) == lengths
        for sequence, index_row in zip(sequential, fai, strict=True):
            text = str(sequence.seq)
            size = len(text)
            width = int(index_row[3])
            intervals = [
                (0, 0),
                (0, 1),
                (0, size),
                (size - 1, size),
                (size, size),
                (size // 2, min(size, size // 2 + 17)),
                (max(0, width - 1), min(size, width + 2)),
            ]
            for _ in range(5):
                start = generator.randrange(size + 1)
                intervals.append((start, generator.randrange(start, size + 1)))
            for start, end in intervals:
                assert indexed.fetch(sequence.id, start, end) == text[start:end], (
                    sequence.id,
                    start,
                    end,
                )
                slice_checks += 1
            assert indexed.fetch(sequence.id) == text, "Full record/case was changed"
            summaries.append(
                {
                    "id": sequence.id,
                    "length": size,
                    "sha256": hashlib.sha256(text.encode("ascii")).hexdigest(),
                }
            )
            normalized = text.upper()
            total_gc += normalized.count("G") + normalized.count("C")
            lowercase += sum(base.islower() for base in text)
            ambiguous += sum(base not in "ACGT" for base in normalized)
    return {
        "sequences": summaries,
        "total_bases": sum(lengths),
        "length_at_least_4": sum(length >= 4 for length in lengths),
        "gc_bases": total_gc,
        "masked_lowercase_bases": lowercase,
        "ambiguous_bases": ambiguous,
        "slice_checks": slice_checks,
    }


def main() -> None:
    """Verify and report the portable scientific result without any BioV service."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    arguments = parser.parse_args()
    _reader_environment()
    prepared = arguments.prepared.resolve()
    record, fasta = _prepared_inventory(prepared)
    summary = _read_sequences(prepared, fasta, record)
    assert _run_readme(prepared) == {
        "records": len(summary["sequences"]),
        "bases": summary["total_bases"],
        "length_at_least_4": summary["length_at_least_4"],
        "literal_GC_fraction": summary["gc_bases"] / summary["total_bases"],
    }
    summary.update(
        {
            "network_blocked": True,
            "inventory_verified": True,
            "readme_example_executed": True,
        }
    )
    print(json.dumps(summary, sort_keys=True))  # noqa: T201 - Standalone acceptance result.


if __name__ == "__main__":
    main()
