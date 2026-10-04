"""Retain full upstream TSV and check its term set, count ratios, and p-values.

Run after ``biov install goatools`` with the same environment configuration.
No downloads, ontology updates, gene-ID inference, or BioV statistics are used.
The output validation record checks selected results; installation and version
evidence must accompany it to establish runtime identity.
"""

import argparse
import csv
import hashlib
import json
import math
import subprocess  # noqa: S404 - fixed upstream CLI acceptance test
import sys
from fractions import Fraction
from pathlib import Path


def fisher_two_sided(a: int, b: int, c: int, d: int) -> Fraction:
    """Enumerate a fixed-margin hypergeometric table using exact integers.

    This independent acceptance-test oracle is not a BioV analysis interface.

    Returns:
        Exact two-sided Fisher probability.
    """
    row = a + b
    successes = a + c
    total = a + b + c + d
    denominator = math.comb(total, row)
    weights = {
        x: math.comb(successes, x) * math.comb(total - successes, row - x)
        for x in range(max(0, row - (total - successes)), min(row, successes) + 1)
    }
    observed = weights[a]
    return Fraction(sum(w for w in weights.values() if w <= observed), denominator)


def require(condition: bool, detail: object) -> None:
    """Raise a validation failure even when Python optimizations are enabled.

    Raises:
        ValueError: If the supplied validation condition fails.
    """
    if not condition:
        raise ValueError(f"GOATOOLS validation failed: {detail}")


def main() -> None:
    """Preserve CLI output, input identities, settings, and independent checks.

    Raises:
        ValueError: If the TSV header is missing or scientific checks fail.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--native-binary", type=Path, help="Use the Rust tools exec route"
    )
    options = parser.parse_args()
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    fixture = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "goatools"
    entry = (
        [str(options.native_binary.resolve()), "tools", "exec", "--no-install"]
        if options.native_binary is not None
        else ["goatools"]
    )
    command = [
        *entry,
        *(["goatools"] if options.native_binary is not None else []),
        "find_enrichment",
        str(fixture / "study.txt"),
        str(fixture / "population.txt"),
        str(fixture / "annotations.id2gos"),
        "--annofmt=id2gos",
        f"--obo={fixture / 'tiny.obo'}",
        "--ns=BP",
        "--alpha=0.05",
        "--method=bonferroni,fdr_bh",
        "--pval=1",
        "--pvalcalc=fisher_scipy_stats",
        f"--outfile={output / 'results.tsv'}",
    ]
    result = subprocess.run(  # noqa: S603 - argv-only, fixed tool and caller output path
        command, capture_output=True, text=True, check=False
    )
    (output / "stdout.txt").write_text(result.stdout)
    (output / "stderr.txt").write_text(result.stderr)
    metadata = {
        "command": command,
        "exitCode": result.returncode,
        "verificationScope": "GO term set, study/population count ratios, selected p-values",
        "runtimeIdentityEvidence": "Requires accompanying installation and version evidence",
        "inputSha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(fixture.iterdir())
            if p.name
            in {"tiny.obo", "study.txt", "population.txt", "annotations.id2gos"}
        },
        "settings": {
            "namespace": "BP",
            "alpha": 0.05,
            "pval": 1,
            "methods": ["bonferroni", "fdr_bh"],
            "pvalcalc": "fisher_scipy_stats",
            "propagate_counts": True,
            "annotation_format": "id2gos",
            "population_n": 10,
            "study_n": 4,
            "tested_terms": 3,
            "hypotheses": "Three BP terms including the propagated root",
            "output_filter": "pval is applied after multiple-testing correction",
        },
    }
    (output / "validation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    result.check_returncode()
    with (output / "results.tsv").open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fieldnames = reader.fieldnames
        if fieldnames is None:
            raise ValueError("GOATOOLS validation failed: Missing TSV header")
        reader.fieldnames = [name.removeprefix("# ") for name in fieldnames]
        records = {}
        for row in reader:
            require(row["GO"] not in records, f"Duplicate GO row: {row['GO']}")
            records[row["GO"]] = row
    require(set(records) == {"GO:0008150", "GO:9000001", "GO:9000002"}, records)
    expected = {
        "GO:0008150": ((4, 4), (10, 10), Fraction(1)),
        "GO:9000001": ((4, 4), (4, 10), fisher_two_sided(4, 0, 0, 6)),
        "GO:9000002": ((0, 4), (6, 10), fisher_two_sided(0, 4, 6, 0)),
    }
    require(expected["GO:9000001"][2] == Fraction(1, 210), expected)
    for go, (study, population, pvalue) in expected.items():
        row = records[go]
        require(row["ratio_in_study"] == f"{study[0]}/{study[1]}", row)
        require(row["ratio_in_pop"] == f"{population[0]}/{population[1]}", row)
        corrections = {
            "p_uncorrected": pvalue,
            "p_bonferroni": min(Fraction(1), 3 * pvalue),
            "p_fdr_bh": Fraction(1) if pvalue == 1 else Fraction(1, 140),
        }
        for field, value in corrections.items():
            # Upstream TSV retains full floating-point values.
            require(math.isclose(float(row[field]), float(value), rel_tol=1e-12), row)
    metadata["verified"] = True
    metadata["exactExpectedLeafPvalues"] = {
        "uncorrected": "1/210",
        "bonferroni": "1/70",
        "fdr_bh": "1/140",
    }
    metadata["resultsSha256"] = hashlib.sha256(
        (output / "results.tsv").read_bytes()
    ).hexdigest()
    (output / "validation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    sys.stdout.write(
        "Verified all three terms, count ratios, Fisher p-values, Bonferroni, and BH\n"
    )


if __name__ == "__main__":
    main()
