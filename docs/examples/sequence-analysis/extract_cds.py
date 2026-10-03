"""Translate explicitly selected, complete standard-code GenBank CDS features."""

import csv
import json
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.SeqFeature import ExactPosition
from Bio.SeqRecord import SeqRecord

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
accessions = parameters.get("accessions")
check(
    "selected_accessions",
    set(parameters) == {"accessions"}
    and isinstance(accessions, list)
    and bool(accessions)
    and all(isinstance(value, str) and value for value in accessions)
    and len(set(accessions)) == len(accessions),
    "Select a nonempty list of unique, versioned GenBank accessions.",
)
try:
    records = SeqIO.to_dict(SeqIO.parse(inputs["genbank"], "genbank"))
except ValueError as exc:
    check("genbank_records", False, str(exc))
    raise
check(
    "genbank_records",
    all(accession in records for accession in accessions),
    f"Requested: {', '.join(accessions)}; available: {', '.join(records)}.",
)

proteins: list[SeqRecord] = []
rows: list[dict[str, str | int]] = []
for accession in accessions:
    record = records[accession]
    features = [feature for feature in record.features if feature.type == "CDS"]
    check(f"{accession}:single_CDS", len(features) == 1, "Require one CDS per record.")
    feature = features[0]
    location = feature.location
    check(
        f"{accession}:complete_location",
        location is not None
        and all(
            isinstance(part.start, ExactPosition)
            and isinstance(part.end, ExactPosition)
            and part.strand in {1, -1}
            and part.ref is None
            for part in location.parts
        ),
        "Require complete local CDS boundaries and an explicit strand.",
    )
    check(
        f"{accession}:standard_code",
        feature.qualifiers.get("transl_table", ["1"]) == ["1"]
        and feature.qualifiers.get("codon_start", ["1"]) == ["1"]
        and "transl_except" not in feature.qualifiers,
        "Use genetic code 1, codon_start 1, and no translation exceptions.",
    )
    check(
        f"{accession}:annotation",
        len(feature.qualifiers.get("protein_id", [])) == 1
        and len(feature.qualifiers.get("translation", [])) == 1,
        "Require one protein_id and one annotated translation.",
    )
    try:
        protein = feature.translate(record.seq, table=1, cds=True)
    except TranslationError as exc:
        check(f"{accession}:translation", False, str(exc))
        raise
    check(
        f"{accession}:translation",
        str(protein) == feature.qualifiers["translation"][0],
        "Computed complete-CDS translation must equal the annotated protein.",
    )
    protein_id = feature.qualifiers["protein_id"][0]
    proteins.append(
        SeqRecord(protein, id=protein_id, description=f"source={accession}")
    )
    rows.append(
        {
            "accession": accession,
            "protein_id": protein_id,
            "cds_location": str(location),
            "genetic_code": 1,
            "length_aa": len(protein),
        }
    )
check(
    "unique_protein_ids",
    len({protein.id for protein in proteins}) == len(proteins),
    "Every selected CDS must retain a distinct protein identifier.",
)
SeqIO.write(proteins, "proteins.fasta", "fasta")
with Path("cds.csv").open("w", newline="") as output:
    writer = csv.DictWriter(output, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
