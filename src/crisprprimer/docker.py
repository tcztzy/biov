"""Run the legacy crisprprimer Docker image and normalize its HTML output."""

import argparse
import json
import subprocess  # noqa: S404 - native software execution
import sys
import tempfile
from pathlib import Path
from typing import cast

import pandas as pd
from pydantic import BaseModel


class CrisprPrimerHit(BaseModel):
    """Validated columns from one legacy HTML report row."""

    rank: int
    id: str
    chrom: str
    pam_start: int
    pam_end: int
    pam_out: str
    spacer_1: str
    spacer_2: str
    score: int


HTML_COLUMN_MAP = {
    "rank": "rank",
    "ID": "id",
    "CHROM": "chrom",
    "PAM start": "pam_start",
    "PAM end": "pam_end",
    "PAM Out": "pam_out",
    "Spacer 1": "spacer_1",
    "Spacer 2": "spacer_2",
    "Score": "score",
}


def find_html_output(output_dir: Path, prefix: str) -> Path:
    """Locate the requested or sole legacy HTML report.

    Returns:
        Path of the native HTML report.

    Raises:
        FileNotFoundError: If no unique report is available.
    """
    expected = output_dir / f"{prefix}.html"
    if expected.exists():
        return expected
    html_files = sorted(output_dir.glob("*.html"))
    if len(html_files) == 1:
        return html_files[0]
    raise FileNotFoundError(f"Could not find a unique HTML report in {output_dir}")


def parse_html_report(report_path: Path) -> list[CrisprPrimerHit]:
    """Read the native HTML table and validate its required columns.

    Returns:
        Validated report rows in original order.

    Raises:
        ValueError: If required columns are absent.
    """
    table = pd.read_html(report_path, header=0)[0].rename(columns=HTML_COLUMN_MAP)
    if "rank" not in table.columns:
        table = table.reset_index(names="rank")
        table["rank"] = table["rank"] + 1
    required_columns = set(CrisprPrimerHit.model_fields)
    missing = required_columns - set(table.columns)
    if missing:
        raise ValueError(f"HTML report is missing columns: {sorted(missing)}")
    return [
        CrisprPrimerHit.model_validate(row)
        for row in cast(pd.DataFrame, table[list(required_columns)]).to_dict(
            orient="records"
        )
    ]


def run_container(
    *,
    reference_path: Path,
    annotation_path: Path,
    region: str,
    prefix: str,
    image: str,
    output_dir: Path,
) -> Path:
    """Run the selected legacy image with explicit mounted inputs.

    Returns:
        Path of the completed report.
    """
    command = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{reference_path.resolve()}:/data/reference.fa:ro",
        "-v",
        f"{annotation_path.resolve()}:/data/annotation.gff:ro",
        "-v",
        f"{output_dir.resolve()}:/app:rw",
        image,
        "-f",
        "/data/reference.fa",
        "-g",
        "/data/annotation.gff",
        "-r",
        region,
        "-h",
        prefix,
    ]
    subprocess.run(  # noqa: S603 - operator-selected programs, literal argv
        command, check=True, capture_output=True, text=True
    )
    return find_html_output(output_dir, prefix)


def build_parser() -> argparse.ArgumentParser:
    """Declare the established command-line arguments.

    Returns:
        The argument parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference", type=Path, required=True, help="Reference FASTA file."
    )
    parser.add_argument(
        "--annotation", type=Path, required=True, help="Annotation GFF/GFF3 file."
    )
    parser.add_argument(
        "--region", required=True, help="Target region in <chrom>:<start>-<end> format."
    )
    parser.add_argument(
        "--prefix", default="A", help="Output prefix passed to crisprprimer."
    )
    parser.add_argument(
        "--image", default="crisprprimer", help="Docker image name to run."
    )
    parser.add_argument(
        "--output-format",
        choices=("json", "csv"),
        default="json",
        help="Normalized output format.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Optional output file. Defaults to stdout.",
    )
    return parser


def main() -> None:
    """Run the command-line workflow with the caller-selected inputs."""
    parser = build_parser()
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="crisprprimer-") as temp_dir:
        report_path = run_container(
            reference_path=args.reference,
            annotation_path=args.annotation,
            region=args.region,
            prefix=args.prefix,
            image=args.image,
            output_dir=Path(temp_dir),
        )
        hits = parse_html_report(report_path)

    if args.output_format == "json":
        payload = json.dumps([hit.model_dump(mode="json") for hit in hits], indent=2)
    else:
        frame = pd.DataFrame([hit.model_dump(mode="json") for hit in hits])
        payload = frame.to_csv(index=False)

    if args.out is None:
        sys.stdout.write(payload + "\n")
    else:
        args.out.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
