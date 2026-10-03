"""Calculate sequence-based average molecular weight and theoretical pI."""

import csv
import json
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

checks: list[dict[str, str | bool]] = []
Path("checks.json").write_text("[]\n")


def check(name: str, passed: bool, detail: str) -> None:
    """Save each data check and stop on a failed precondition.

    Raises:
        ValueError: When the scientific input fails the stated check.
    """
    checks.append({"name": name, "passed": passed, "detail": detail})
    Path("checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    if not passed:
        raise ValueError(f"{name}: {detail}")


inputs = json.loads(Path(sys.argv[1]).read_text())
parameters = json.loads(Path(sys.argv[2]).read_text())
check("fixed_method", parameters == {}, "This example accepts no method parameters.")
try:
    proteins = list(SeqIO.parse(inputs["proteins"], "fasta"))
except ValueError as exc:
    check("protein_records", False, str(exc))
    raise
check("protein_records", bool(proteins), "Require at least one FASTA protein record.")
check(
    "unique_protein_ids",
    len({protein.id for protein in proteins}) == len(proteins),
    "Protein identifiers must be unique.",
)
rows: list[dict[str, str | int | float]] = []
for record in proteins:
    sequence = str(record.seq)
    check(
        f"{record.id}:canonical_protein",
        bool(sequence) and set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY"),
        "Require nonempty canonical amino acids; stops and ambiguous residues fail.",
    )
    analysis = ProteinAnalysis(sequence, monoisotopic=False)
    rows.append(
        {
            "protein_id": record.id,
            "length_aa": len(sequence),
            "molecular_weight_da": analysis.molecular_weight(),
            "theoretical_pi": analysis.isoelectric_point(),
        }
    )
with Path("properties.csv").open("w", newline="") as output:
    writer = csv.DictWriter(output, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
