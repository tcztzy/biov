"""Inspect the example GCF_000005845.2 package without BioV or network calls.

This is example-specific, not a general GFF validator. Its first CDS is unspliced
and its coordinates do not cross a circular origin. Counts describe this snapshot.
"""

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any


def require(condition: bool, message: str) -> None:
    """Fail a verification independently of Python optimization settings.

    Raises:
        ValueError: If the supplied verification condition is false.
    """
    if not condition:
        raise ValueError(message)


def safe_file(base: Path, name: str) -> Path:
    """Resolve a native relative filename without traversal or escaping links.

    Returns:
        An existing file confined to the supplied trusted package directory.
    """
    require(
        bool(name)
        and "\\" not in name
        and not name.startswith("/")
        and ":" not in name.split("/", maxsplit=1)[0]
        and all(part not in {"", ".", ".."} for part in name.split("/")),
        f"Unsafe relative path: {name!r}",
    )
    target = base / name
    require(
        target.is_file() and target.resolve().is_relative_to(base.resolve()),
        f"Unsafe or missing file: {name!r}",
    )
    return target


def read_fasta(path: Path) -> dict[str, str]:
    """Read small FASTA fixtures and reject duplicate sequence identifiers.

    Returns:
        Exact sequence strings keyed by the first whitespace-delimited header ID.
    """
    records: dict[str, str] = {}
    name: str | None = None
    sequence: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(sequence)
            fields = line[1:].split()
            require(bool(fields), "Empty FASTA header")
            name = fields[0]
            require(name not in records, f"Duplicate FASTA ID: {name}")
            sequence = []
        elif line:
            require(name is not None, "Sequence precedes first FASTA header")
            sequence.append(line)
    if name is not None:
        records[name] = "".join(sequence)
    require(bool(records), f"Empty FASTA file: {path.name}")
    return records


def inspect(root: Path) -> dict[str, Any]:
    """Verify and inspect this E. coli acquisition using its native catalog.

    Returns:
        Check results and a small complete sequence example from this acquisition.

    Raises:
        ValueError: If required package checks or example-specific assumptions fail.
    """
    root = root.resolve()
    data = root / "ncbi_dataset/data"
    require(data.resolve().is_relative_to(root), "Data directory escapes package")
    catalog = json.loads(
        safe_file(root, "ncbi_dataset/data/dataset_catalog.json").read_text()
    )
    assemblies = [
        item
        for item in catalog["assemblies"]
        if item.get("accession") == "GCF_000005845.2"
    ]
    require(len(assemblies) == 1, "Expected exactly one example assembly")
    paths: dict[str, Path] = {}
    for entry in assemblies[0]["files"]:
        kind = entry["fileType"]
        require(kind not in paths, f"Duplicate example file type: {kind}")
        paths[kind] = safe_file(data, entry["filePath"])
    seen_paths: set[str] = set()
    for group in catalog["assemblies"]:
        for entry in group["files"]:
            name = entry["filePath"]
            require(name not in seen_paths, f"Duplicate catalog path: {name}")
            seen_paths.add(name)
            require(
                safe_file(data, name).stat().st_size
                == int(entry["uncompressedLengthBytes"]),
                f"Catalog size mismatch: {name}",
            )
    checksum_paths: set[str] = set()
    for line in safe_file(root, "md5sum.txt").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        require(name not in checksum_paths, f"Duplicate checksum path: {name}")
        checksum_paths.add(name)
        actual = hashlib.md5(safe_file(root, name).read_bytes(), usedforsecurity=False)
        require(actual.hexdigest() == digest, f"MD5 mismatch: {name}")
    expected_checksums = {f"ncbi_dataset/data/{name}" for name in seen_paths} | {
        "ncbi_dataset/data/dataset_catalog.json"
    }
    require(
        checksum_paths == expected_checksums,
        "Example checksums must cover the native catalog and every catalog file",
    )

    genome = read_fasta(paths["GENOMIC_NUCLEOTIDE_FASTA"])
    proteins = read_fasta(paths["PROTEIN_FASTA"])
    cds = read_fasta(paths["CDS_NUCLEOTIDE_FASTA"])
    features: Counter[str] = Counter()
    gff_proteins: set[str] = set()
    first: tuple[list[str], dict[str, str]] | None = None
    for line in paths["GFF3"].read_text().splitlines():
        if line.startswith("#") or not line:
            continue
        fields = line.split("\t")
        require(len(fields) == 9, "Expected nine GFF3 columns")
        features[fields[2]] += 1
        require(fields[0] in genome, "GFF sequence ID missing from genome")
        require(
            1 <= int(fields[3]) <= int(fields[4]) <= len(genome[fields[0]]),
            "Example GFF coordinates outside genome",
        )
        attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
        if "protein_id" in attrs:
            gff_proteins.add(attrs["protein_id"])
        if first is None and fields[2] == "CDS":
            first = fields, attrs
    require(gff_proteins == set(proteins), "GFF protein IDs differ from protein FASTA")
    report = json.loads(paths["SEQUENCE_REPORT"].read_text().splitlines()[0])
    require(
        report["refseqAccession"] in genome, "Sequence report ID missing from genome"
    )
    require(
        report["length"] == len(genome[report["refseqAccession"]]),
        "Sequence length mismatch",
    )
    if first is None:
        raise ValueError("Example has no CDS")
    fields, attrs = first
    require(
        fields[0] == "NC_000913.3"
        and fields[3:5] == ["190", "255"]
        and fields[6:8] == ["+", "0"]
        and attrs.get("gene") == "thrL"
        and attrs.get("protein_id") == "NP_414542.1",
        "First CDS is not the documented plus-strand thrL example",
    )
    sequence = genome[fields[0]][int(fields[3]) - 1 : int(fields[4])]
    require(
        set(sequence) <= set("ACGT"),
        "Example CDS contains unsupported ambiguity symbols",
    )
    require(
        cds.get("lcl|NC_000913.3_cds_NP_414542.1_1") == sequence,
        "First example CDS differs from genome",
    )
    return {
        "all_native_md5_pass": True,
        "native_md5_entries": len(checksum_paths),
        "catalog_sizes_match": True,
        "fasta_records": {
            "genome": len(genome),
            "protein": len(proteins),
            "cds": len(cds),
        },
        "genome_lengths": {key: len(value) for key, value in genome.items()},
        "gff_feature_counts": dict(features),
        "gff_seqids_match_genome": True,
        "gff_protein_ids_exactly_match_protein_fasta": True,
        "sequence_report": report,
        "first_cds": {
            "protein_id": attrs["protein_id"],
            "gene": attrs.get("gene"),
            "start_1based": int(fields[3]),
            "end_inclusive": int(fields[4]),
            "strand": fields[6],
            "nucleotides": sequence,
            "matches_cds_fasta": True,
            "protein_sequence": proteins[attrs["protein_id"]],
        },
        "rna_fasta_present": "RNA_NUCLEOTIDE_FASTA" in paths,
    }


def main() -> None:
    """Run the explicit example package verifier."""
    parser = argparse.ArgumentParser(
        description="Offline stdlib inspection of the GCF_000005845.2 example; not a general GFF validator."
    )
    parser.add_argument(
        "package_root",
        type=Path,
        help="Extracted NCBI package root containing README.md",
    )
    args = parser.parse_args()
    sys.stdout.write(json.dumps(inspect(args.package_root), indent=2) + "\n")


if __name__ == "__main__":
    main()
