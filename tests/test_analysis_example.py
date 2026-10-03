"""Check the real GenBank-to-protein example independently of MCP execution."""

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "docs/examples/sequence-analysis"
SOURCE = ROOT / "tests/data/cor6_6.gb"
ACCESSIONS = ["X55053.1", "X62281.1", "M81224.1", "L31939.1", "AF297471.1"]
PROTEINS = ["CAA38894.1", "CAA44171.1", "AAA32993.1", "AAA91051.1", "AAG13407.1"]


def run_example(
    directory: Path,
    script: str,
    inputs: dict[str, Path],
    parameters: dict[str, object],
) -> subprocess.CompletedProcess[str]:
    """Run one example script with the analysis worker's file-based contract.

    Returns:
        Captured process status and diagnostics.
    """
    directory.mkdir()
    input_path = directory / "inputs.json"
    parameter_path = directory / "parameters.json"
    input_path.write_text(
        json.dumps({name: str(path) for name, path in inputs.items()})
    )
    parameter_path.write_text(json.dumps(parameters))
    return subprocess.run(
        [sys.executable, str(EXAMPLE / script), str(input_path), str(parameter_path)],
        cwd=directory,
        capture_output=True,
        text=True,
        check=False,
    )


def test_real_complete_cds_translation_and_protein_properties(tmp_path: Path) -> None:
    """Retain all five records, joined CDS coordinates, and fixed scientific values."""
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == (
        "01b4e193b71344752a96e2b37e886118a33c07b16c3684741c31eb3efeaf60b8"
    )
    extraction = tmp_path / "extraction"
    result = run_example(
        extraction, "extract_cds.py", {"genbank": SOURCE}, {"accessions": ACCESSIONS}
    )
    assert result.returncode == 0, result.stderr
    proteins = list(SeqIO.parse(extraction / "proteins.fasta", "fasta"))
    assert [protein.id for protein in proteins] == PROTEINS
    assert [len(protein) for protein in proteins] == [66, 67, 65, 65, 65]
    with (extraction / "cds.csv").open() as source:
        cds = list(csv.DictReader(source))
    assert [row["accession"] for row in cds] == ACCESSIONS
    assert cds[0]["cds_location"] == "[49:250](+)"
    assert cds[1]["cds_location"] == ("join{[103:160](+), [319:390](+), [503:579](+)}")
    assert cds[-1]["cds_location"] == "join{[0:54](+), [240:309](+), [422:497](+)}"
    assert all(row["genetic_code"] == "1" for row in cds)

    properties = tmp_path / "properties"
    result = run_example(
        properties,
        "protein_properties.py",
        {"proteins": extraction / "proteins.fasta"},
        {},
    )
    assert result.returncode == 0, result.stderr
    with (properties / "properties.csv").open() as source:
        rows = list(csv.DictReader(source))
    assert [row["protein_id"] for row in rows] == PROTEINS
    assert [float(row["molecular_weight_da"]) for row in rows] == pytest.approx(
        [6551.1852, 7407.5077, 6552.2611, 6604.4217, 6536.2617]
    )
    assert [float(row["theoretical_pi"]) for row in rows] == pytest.approx(
        [9.1007600784, 10.9675710678, 9.1576856613, 9.1576856613, 9.1576856613]
    )
    for directory in (extraction, properties):
        checks = json.loads((directory / "checks.json").read_text())
        assert checks and all(check["passed"] for check in checks)


@pytest.mark.parametrize("failure", ["partial", "annotation", "genetic_code"])
def test_cds_scientific_failures_are_recorded(tmp_path: Path, failure: str) -> None:
    """Reject incomplete CDSs, mismatched translation, and unsupported genetic codes."""
    source = tmp_path / "input.gb"
    content = SOURCE.read_text()
    accessions = ["X55053.1"]
    expected = "translation"
    if failure == "partial":
        accessions = ["AJ237582.1"]
        expected = "complete_location"
    elif failure == "annotation":
        content = content.replace('/translation="MSET', '/translation="ASET', 1)
    else:
        content = content.replace(
            "/codon_start=1", "/codon_start=1\n                     /transl_table=2", 1
        )
        expected = "standard_code"
    source.write_text(content)
    directory = tmp_path / "run"
    result = run_example(
        directory, "extract_cds.py", {"genbank": source}, {"accessions": accessions}
    )
    assert result.returncode != 0
    checks = json.loads((directory / "checks.json").read_text())
    assert checks[-1]["passed"] is False
    assert checks[-1]["name"].endswith(expected)
    assert not (directory / "proteins.fasta").exists()


def test_protein_properties_rejects_ambiguous_residues(tmp_path: Path) -> None:
    """Do not silently remove ambiguous residues or stop symbols before calculation."""
    source = tmp_path / "input.fasta"
    source.write_text(">invalid\nMX*\n")
    directory = tmp_path / "run"
    result = run_example(directory, "protein_properties.py", {"proteins": source}, {})
    assert result.returncode != 0
    checks = json.loads((directory / "checks.json").read_text())
    assert checks[-1]["name"] == "invalid:canonical_protein"
    assert checks[-1]["passed"] is False
    assert not (directory / "properties.csv").exists()
