"""Independent acceptance consumer copied out of the checkout before execution.

Only the Python standard library is used for inventories, synthetic data and the
RefSeq example. The real PDB check additionally uses unmodified Biopython. No
BioV module, running service, catalog database or original import path is used.
This helper is not a BioV runtime API or a general biological-format validator.
"""

import argparse
import contextlib
import errno
import gzip
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import os
import re
import runpy
import socket
import struct
import sys
from pathlib import Path
from typing import Any


def _reader_environment(mode: str) -> None:
    """Prove interpreter isolation and fail if kernel networking remains usable.

    Raises:
        AssertionError: If socket creation succeeds despite the required filter.
    """
    assert os.environ["PATH"] == os.environ["PYTHONPATH"] == ""
    assert sys.prefix != sys.base_prefix, "Reader must use a separate virtualenv"
    assert importlib.util.find_spec("biov") is None, "BioV must not be installed"
    packages = {
        distribution.metadata["Name"].lower()
        for distribution in importlib.metadata.distributions()
    }
    expected = {"biopython", "numpy"} if mode == "pdb" else {"pyarrow"}
    assert packages == expected, (packages, expected)
    # This call reaches the OS before the additional Python audit hook is added.
    # Merely omitting an HTTP call, clearing proxies or monkeypatching requests
    # would not establish that either Rust or native readers cannot use network.
    for family, kind in (
        (socket.AF_INET, socket.SOCK_STREAM),
        (socket.AF_INET6, socket.SOCK_DGRAM),
        (socket.AF_UNIX, socket.SOCK_STREAM),
    ):
        try:
            sock = socket.socket(family, kind)
        except PermissionError as error:
            assert error.errno == errno.EPERM
        else:
            sock.close()
            raise AssertionError("Kernel network filter is not active")

    def offline_only(event: str, args: tuple[Any, ...]) -> None:
        """Reject accidental dynamic BioV/database dependencies and subprocesses.

        Raises:
            AssertionError: If a forbidden dependency or operation is attempted.
        """
        if event.startswith(("socket.", "subprocess.", "sqlite3.")) or event in {
            "os.system",
            "os.exec",
            "os.posix_spawn",
            "os.fork",
            "os.forkpty",
        }:
            raise AssertionError(f"Reader must remain independent/offline: {event}")
        if event == "import" and args[0].split(".")[0] == "biov":
            raise AssertionError(f"Forbidden reader dependency: {args[0]}")

    sys.addaudithook(offline_only)


def _safe_file(root: Path, name: str) -> Path:
    """Return a confined ordinary file from a portable relative receipt path."""
    assert name and "\\" not in name and not Path(name).is_absolute(), name
    assert all(part not in {"", ".", ".."} for part in name.split("/")), name
    path = root / name
    assert path.is_file() and not path.is_symlink(), name
    assert path.resolve().is_relative_to(root.resolve()), name
    return path


def _check_inventory(snapshot: Path) -> dict[str, Any]:
    """Verify the complete tree, documented canonical digest and SHA-256 list.

    Returns:
        A plain-JSON receipt whose identities agree with every saved source byte.
    """
    receipt = json.loads(_safe_file(snapshot, "acquisition.json").read_text())
    assert receipt["schema_version"] == 1
    assert receipt["origin"] == "local_copy"
    assert receipt["registered_at_unix_seconds"] > 0
    for field in ("acquired_at", "source_url", "acquisition_tool"):
        assert receipt[field] is None, f"Unknown acquisition fact was invented: {field}"
    source = snapshot / "source"
    assert source.is_dir() and not source.is_symlink()
    entries = receipt["inventory"]
    names = [entry["path"] for entry in entries]
    assert names == sorted(set(names), key=lambda name: name.encode("utf-8"))
    actual_names = {path.relative_to(source).as_posix() for path in source.rglob("*")}
    assert set(names) == actual_names, "Inventory must cover all files and directories"
    digest = hashlib.sha256(b"biov-native-tree-v1\0")
    checksums = {}
    for entry in entries:
        name = entry["path"]
        raw = name.encode("utf-8")
        assert not (source / name).is_symlink(), name
        if entry["kind"] == "directory":
            assert (source / name).is_dir()
            assert entry["bytes"] == 0 and entry["sha256"] is None
            kind, size, content_digest = b"D", 0, bytes(32)
        else:
            assert entry["kind"] == "file", entry
            path = _safe_file(source, name)
            with path.open("rb") as stream:
                sha256 = hashlib.file_digest(stream, "sha256").hexdigest()
            size = path.stat().st_size
            assert entry["bytes"] == size and entry["sha256"] == sha256, name
            kind, content_digest = b"F", bytes.fromhex(sha256)
            checksums["source/" + name] = sha256
        digest.update(kind + struct.pack(">I", len(raw)) + raw)
        digest.update(struct.pack(">Q", size) + content_digest)
    assert receipt["source_content_sha256"] == digest.hexdigest()
    assert receipt["snapshot_id"] == "sha256-" + digest.hexdigest()
    assert snapshot.name == receipt["snapshot_id"]
    lines = _safe_file(snapshot, "checksums.sha256").read_text().splitlines()
    observed = {}
    for line in lines:
        expected, name = line.split(maxsplit=1)
        assert name not in observed
        observed[name] = expected
    assert observed == checksums
    readme_bytes = _safe_file(snapshot, "README.md").read_bytes()
    assert hashlib.sha256(readme_bytes).hexdigest() == receipt["readme_sha256"]
    readme = readme_bytes.decode("utf-8")
    assert "acquisition.json" in readme and "checksums.sha256" in readme
    for representation, paths in receipt["representations"].items():
        assert representation in readme, (
            "README must list every analysis representation"
        )
        for name in paths:
            _safe_file(source, name)
            assert "source/" + name in readme, "README must give direct ordinary paths"
    for name in receipt["native_metadata"]:
        _safe_file(source, name)
    return receipt


def _run_readme(snapshot: Path) -> str:
    """Return output from the emitted, independently executable README example."""
    readme = _safe_file(snapshot, "README.md").read_text()
    examples = re.findall(
        r"^```python\s*\n(.*?)^```\s*$", readme, re.MULTILINE | re.DOTALL
    )
    assert len(examples) == 1, "README must provide a runnable native-reader example"
    previous = Path.cwd()
    output = io.StringIO()
    try:
        os.chdir(snapshot)
        with contextlib.redirect_stdout(output):
            exec(  # noqa: S102 - Run the trusted emitted README, not provider code.
                compile(examples[0], "snapshot README example", "exec"),
                {"__name__": "__main__"},
            )
    finally:
        os.chdir(previous)
    return output.getvalue().strip()


def _inspect_synthetic(source: Path) -> dict[str, Any]:
    """Read every synthetic FASTA record and calculate a direct sequence summary.

    Returns:
        Full independently calculated counts and sequence lengths.
    """
    catalog = json.loads(
        (source / "ncbi_dataset/data/dataset_catalog.json").read_text()
    )
    group = next(group for group in catalog["assemblies"] if group.get("accession"))
    entry = next(
        entry
        for entry in group["files"]
        if entry["fileType"] == "GENOMIC_NUCLEOTIDE_FASTA"
    )
    path = _safe_file(source / "ncbi_dataset/data", entry["filePath"])
    sequences = {}
    identifier = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            identifier = line[1:].split()[0]
            assert identifier not in sequences
            sequences[identifier] = ""
        elif line:
            assert identifier is not None
            sequences[identifier] += line
    assert sequences == {"synthetic_a": "ACGTNN", "synthetic_b": "GGCC"}
    assert {name: seq for name, seq in sequences.items() if "N" not in seq} == {
        "synthetic_b": "GGCC"
    }
    return {
        "fasta_records": len(sequences),
        "total_bases": sum(map(len, sequences.values())),
        "unambiguous_records": 1,
    }


def _inspect_refseq(source: Path, reader: Path) -> dict[str, Any]:
    """Use the existing independent example reader on all real RefSeq records.

    Returns:
        Verified full-record counts and the native thrL coordinate cross-check.
    """
    summary = runpy.run_path(str(reader))["inspect"](source)
    assert summary["all_native_md5_pass"] is True
    assert summary["native_md5_entries"] == 7
    assert summary["catalog_sizes_match"] is True
    assert summary["fasta_records"] == {"genome": 1, "protein": 4300, "cds": 4318}
    assert summary["genome_lengths"] == {"NC_000913.3": 4641652}
    assert summary["gff_feature_counts"] == {
        "region": 1,
        "gene": 4506,
        "CDS": 4340,
        "mobile_genetic_element": 50,
        "ncRNA": 108,
        "exon": 216,
        "rRNA": 22,
        "tRNA": 86,
        "pseudogene": 145,
        "sequence_feature": 48,
        "origin_of_replication": 1,
    }
    assert summary["gff_seqids_match_genome"] is True
    assert summary["gff_protein_ids_exactly_match_protein_fasta"] is True
    assert summary["rna_fasta_present"] is False
    assert summary["first_cds"]["gene"] == "thrL"
    assert summary["first_cds"]["protein_sequence"] == "MKRISTTITTTITITTGNGAG"
    assert summary["first_cds"]["matches_cds_fasta"] is True
    assert (source / "README.md").is_file()
    assert (source / "ncbi_dataset/data/assembly_data_report.jsonl").is_file()
    return summary


def _inspect_pdb(source: Path) -> dict[str, Any]:
    """Read the complete real mmCIF with stock Biopython and native field names.

    Returns:
        Independently counted models, residues and atoms with native metadata.
    """
    import Bio
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from Bio.PDB.MMCIFParser import MMCIFParser

    assert Bio.__version__ == "1.88"
    assert {path.name for path in source.iterdir()} == {
        "1crn.cif",
        "1crn.cif.gz",
        "1crn.entry.json",
    }, "Do not register investigation scripts or derived inspection reports"
    cif = _safe_file(source, "1crn.cif")
    assert gzip.decompress((source / "1crn.cif.gz").read_bytes()) == cif.read_bytes()
    structure = MMCIFParser(QUIET=True).get_structure("1CRN", str(cif))
    assert structure is not None
    metadata = MMCIF2Dict(str(cif))
    summary = {
        "entry": metadata["_entry.id"][0],
        "models": len(structure),
        "chains": [chain.id for chain in structure[0]],
        "residues": sum(1 for _ in structure.get_residues()),
        "atoms": sum(1 for _ in structure.get_atoms()),
        "method": metadata["_exptl.method"][0],
        "resolution_angstrom": float(metadata["_refine.ls_d_res_high"][0]),
        "revision": ".".join(
            metadata[f"_pdbx_audit_revision_history.{version}_revision"][-1]
            for version in ("major", "minor")
        ),
        "revision_date": metadata["_pdbx_audit_revision_history.revision_date"][-1],
    }
    assert summary == {
        "entry": "1CRN",
        "models": 1,
        "chains": ["A"],
        "residues": 46,
        "atoms": 327,
        "method": "X-RAY DIFFRACTION",
        "resolution_angstrom": 1.5,
        "revision": "1.5",
        "revision_date": "2024-10-30",
    }
    entry = json.loads((source / "1crn.entry.json").read_text())
    assert entry["rcsb_id"] == summary["entry"]
    assert entry["rcsb_entry_info"]["deposited_atom_count"] == summary["atoms"]
    return summary


def main() -> None:
    """Run a complete independent analysis of one copied native snapshot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("synthetic", "refseq", "pdb"))
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--refseq-reader", type=Path)
    args = parser.parse_args()
    _reader_environment(args.mode)
    receipt = _check_inventory(args.snapshot)
    source = args.snapshot / "source"
    example_output = _run_readme(args.snapshot)
    assert (
        example_output
        == {"synthetic": "2 10", "refseq": "1 4641652", "pdb": "46 327"}[args.mode]
    ), example_output
    if args.mode == "synthetic":
        summary = _inspect_synthetic(source)
    elif args.mode == "refseq":
        assert args.refseq_reader is not None
        summary = _inspect_refseq(source, args.refseq_reader)
        assert "rna_fasta" in receipt["unavailable_representations"]
    else:
        summary = _inspect_pdb(source)
        assert receipt["declaration"]["provider"] == "pdb"
        assert receipt["scope"] == "entry"
        assert receipt["validation"]["method"] == (
            "caller_declared_pdb_identity_scope_and_representations"
        )
    assert not any(name == "biov" or name.startswith("biov.") for name in sys.modules)
    # Stock Bio.File imports sqlite3 for optional indexing; the audit hook
    # forbids opening any database, and mmCIF parsing does not need one.
    if args.mode != "pdb":
        assert "sqlite3" not in sys.modules
    sys.stdout.write(
        json.dumps(
            {"inventory_verified": True, "network_blocked": True, "summary": summary}
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
