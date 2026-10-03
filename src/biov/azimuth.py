"""Score Azimuth 30mers with a Docker-first backend and a pinned fork fallback."""

import argparse
import csv
import json
import os
import shutil
import subprocess  # noqa: S404 - native software execution
import sys
import tempfile
import textwrap
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ValidationError, field_validator, model_validator

from .config import settings

DEFAULT_DOCKER_IMAGE = "azimuth:latest"
DEFAULT_FORK_URL = "https://github.com/milescsmith/Azimuth.git"
DEFAULT_FORK_REF = "3fa94432e924d51e4741efb5d91facb582b4fbfe"
DEFAULT_CACHE_ROOT = settings.home / "azimuth"
FORK_PYTHON = "3.10"
FORK_RUNTIME_DEPENDENCIES = (
    "numpy<2",
    "pandas<3",
    "scikit-learn==1.3.2",
    "scipy<1.12",
    "biopython==1.81",
    "dill==0.3.7",
    "typer",
    "rich",
    "loguru",
    "better-exceptions",
    "tqdm",
)
DOCKER_RUNNER = textwrap.dedent(
    """\
    import json
    import sys

    import numpy as np
    from azimuth.model_comparison import predict

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    records = json.loads(open(input_path).read())
    scores = [None] * len(records)

    with_positions = [
        index
        for index, record in enumerate(records)
        if record["aa_cut"] is not None and record["percent_peptide"] is not None
    ]
    without_positions = [
        index
        for index, record in enumerate(records)
        if record["aa_cut"] is None and record["percent_peptide"] is None
    ]

    if with_positions:
        seq = np.array([records[index]["sequence"] for index in with_positions])
        aa_cut = np.array([records[index]["aa_cut"] for index in with_positions], dtype=float)
        percent_peptide = np.array(
            [records[index]["percent_peptide"] for index in with_positions], dtype=float
        )
        predicted = predict(seq, aa_cut=aa_cut, percent_peptide=percent_peptide)
        for index, score in zip(with_positions, predicted):
            scores[index] = float(score)

    if without_positions:
        seq = np.array([records[index]["sequence"] for index in without_positions])
        predicted = predict(seq, aa_cut=None, percent_peptide=None)
        for index, score in zip(without_positions, predicted):
            scores[index] = float(score)

    with open(output_path, "w") as handle:
        json.dump(scores, handle)
    """
)
FORK_RUNNER = textwrap.dedent(
    """\
    import json
    import sys
    from pathlib import Path

    import numpy as np
    import pandas as pd
    from dill import loads

    source_dir = Path(sys.argv[1])
    input_path = Path(sys.argv[2])
    output_path = Path(sys.argv[3])
    sys.path.insert(0, str(source_dir / "src"))

    from azimuth.features.featurization import featurize_data
    from azimuth.util import concatenate_feature_sets

    def main():
        records = json.loads(input_path.read_text(encoding="utf-8"))
        scores = [None] * len(records)
        model_root = source_dir / "src" / "azimuth" / "azure_models"

        def score_batch(indexes, model_name, use_position):
            try:
                model, learn_options = loads((model_root / model_name).read_bytes())
            except Exception as exc:
                raise RuntimeError(
                    "The pinned milescsmith/Azimuth fallback could not load its bundled model "
                    f"{model_name}. The archived fork ships version-sensitive pickles; restore the "
                    "Docker backend for canonical scoring on this host."
                ) from exc

            learn_options["V"] = 2
            seq = np.array([records[index]["sequence"] for index in indexes])
            x_df = pd.DataFrame(
                columns=["30mer", "Strand"],
                data=list(zip(seq, ["NA" for _ in range(len(seq))], strict=True)),
            )
            if use_position:
                gene_position = pd.DataFrame(
                    columns=["Percent Peptide", "Amino Acid Cut position"],
                    data=list(
                        zip(
                            [records[index]["percent_peptide"] for index in indexes],
                            [records[index]["aa_cut"] for index in indexes],
                            strict=True,
                        )
                    ),
                )
            else:
                gene_position = pd.DataFrame(
                    columns=["Percent Peptide", "Amino Acid Cut position"],
                    data=list(
                        zip(
                            np.ones(seq.shape[0]) * -1,
                            np.ones(seq.shape[0]) * -1,
                            strict=True,
                        )
                    ),
                )

            feature_sets = featurize_data(
                x_df,
                learn_options,
                pd.DataFrame(),
                gene_position,
                pam_audit=True,
                length_audit=False,
            )
            inputs, *_ = concatenate_feature_sets(feature_sets)
            predicted = model.predict(inputs)
            for index, score in zip(indexes, predicted):
                scores[index] = float(score)

        with_positions = [
            index
            for index, record in enumerate(records)
            if record["aa_cut"] is not None and record["percent_peptide"] is not None
        ]
        without_positions = [
            index
            for index, record in enumerate(records)
            if record["aa_cut"] is None and record["percent_peptide"] is None
        ]

        if with_positions:
            score_batch(with_positions, "V3_model_full.pickle", True)
        if without_positions:
            score_batch(without_positions, "V3_model_nopos.pickle", False)

        output_path.write_text(json.dumps(scores), encoding="utf-8")

    if __name__ == "__main__":
        try:
            main()
        except RuntimeError as exc:
            raise SystemExit(str(exc))
    """
)


class SequenceRecord(BaseModel):
    """Validated Azimuth context and optional coding-position metadata."""

    record_id: str | None = None
    sequence: str
    aa_cut: float | None = None
    percent_peptide: float | None = None

    @field_validator("record_id", mode="before")
    @classmethod
    def normalize_record_id(cls, value: Any) -> str | None:
        """Normalize an optional caller label.

        Returns:
            Normalized label or None.
        """
        if value in (None, ""):
            return None
        return str(value)

    @field_validator("sequence", mode="before")
    @classmethod
    def normalize_sequence(cls, value: Any) -> str:
        """Normalize and validate an unambiguous 30-base context.

        Returns:
            Uppercase context of exactly 30 unambiguous bases.

        Raises:
            TypeError: If the input is not text.
            ValueError: If length or alphabet is unsupported.
        """
        if not isinstance(value, str):
            raise TypeError("sequence must be a string")
        sequence = value.strip().upper()
        if len(sequence) != 30:
            raise ValueError("sequence must be exactly 30 nt long")
        if set(sequence) - {"A", "C", "G", "T"}:
            raise ValueError("sequence must contain only A/C/G/T")
        return sequence

    @field_validator("aa_cut", "percent_peptide", mode="before")
    @classmethod
    def normalize_optional_float(cls, value: Any) -> float | None:
        """Parse an optional numeric input field.

        Returns:
            Parsed float or None.
        """
        if value in (None, ""):
            return None
        return float(value)

    @model_validator(mode="after")
    def require_complete_position_pair(self) -> "SequenceRecord":
        """Require both coding-position fields together.

        Returns:
            The validated record.

        Raises:
            ValueError: If only one position field is supplied.
        """
        has_aa_cut = self.aa_cut is not None
        has_percent_peptide = self.percent_peptide is not None
        if has_aa_cut != has_percent_peptide:
            raise ValueError("aa_cut and percent_peptide must be provided together")
        return self


class ScoredRecord(SequenceRecord):
    """One input with its computed score and actual backend."""

    azimuth_score: float
    backend: Literal["docker", "fork"]


def docker_image_available(image: str) -> bool:
    """Check whether the operator has the requested local Docker image.

    Returns:
        Whether Docker can inspect that image locally.
    """
    if shutil.which("docker") is None:
        return False
    command = ["docker", "image", "inspect", image]
    result = subprocess.run(  # noqa: S603 - operator-selected programs, literal argv
        command, check=False, capture_output=True, text=True
    )
    return result.returncode == 0


def choose_backend(
    preferred: Literal["auto", "docker", "fork"], image: str
) -> Literal["docker", "fork"]:
    """Honor an explicit backend or select the available local image.

    Returns:
        The selected docker or fork backend name.
    """
    if preferred == "docker":
        return "docker"
    if preferred == "fork":
        return "fork"
    return "docker" if docker_image_available(image) else "fork"


def normalize_input_mapping(raw: dict[str, Any]) -> dict[str, Any]:
    """Map the documented input column aliases.

    Returns:
        Fields accepted by SequenceRecord.
    """
    return {
        "record_id": raw.get("record_id", raw.get("id")),
        "sequence": raw.get("sequence", raw.get("30mer")),
        "aa_cut": raw.get("aa_cut"),
        "percent_peptide": raw.get("percent_peptide"),
    }


def load_records_from_path(path: Path) -> list[SequenceRecord]:
    """Read and validate a JSON, JSONL, CSV, or TSV batch.

    Returns:
        Validated records in file order.

    Raises:
        ValueError: If the input extension is unsupported.
    """
    suffix = path.suffix.lower()
    if suffix == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            payload = [payload]
    elif suffix == ".jsonl":
        payload = [
            json.loads(line)
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    elif suffix in {".csv", ".tsv"}:
        delimiter = "," if suffix == ".csv" else "\t"
        with path.open(encoding="utf-8", newline="") as handle:
            payload = list(csv.DictReader(handle, delimiter=delimiter))
    else:
        raise ValueError(f"Unsupported input file format: {path.suffix}")
    return [
        SequenceRecord.model_validate(normalize_input_mapping(item)) for item in payload
    ]


def load_records_from_args(args: argparse.Namespace) -> list[SequenceRecord]:
    """Validate the selected file or command-line contexts.

    Returns:
        Validated records in argument order.

    Raises:
        ValueError: If no input or an ambiguous position batch is supplied.
    """
    if args.input is not None:
        return load_records_from_path(args.input)
    if not args.sequence:
        raise ValueError("Provide either --input or at least one --sequence")
    if args.aa_cut is None and args.percent_peptide is None:
        return [SequenceRecord(sequence=sequence) for sequence in args.sequence]
    if len(args.sequence) != 1:
        raise ValueError(
            "Positional scoring through CLI flags supports exactly one --sequence"
        )
    return [
        SequenceRecord(
            sequence=args.sequence[0],
            aa_cut=args.aa_cut,
            percent_peptide=args.percent_peptide,
        )
    ]


def write_helper_script(path: Path, contents: str) -> None:
    """Save the exact backend runner in its temporary directory."""
    path.write_text(contents, encoding="utf-8")


def write_record_batch(path: Path, records: list[SequenceRecord]) -> None:
    """Serialize validated inputs for the selected scoring process."""
    path.write_text(
        json.dumps([record.model_dump(mode="json") for record in records], indent=2),
        encoding="utf-8",
    )


def build_docker_command(
    *,
    image: str,
    work_dir: Path,
) -> list[str]:
    """Build the existing container invocation with literal arguments.

    Returns:
        Literal Docker argv.
    """
    return [
        "docker",
        "run",
        "--rm",
        "--platform",
        "linux/amd64",
        "-v",
        f"{work_dir.resolve()}:/work",
        image,
        "/work/runner.py",
        "/work/input.json",
        "/work/output.json",
    ]


def run_docker_backend(records: list[SequenceRecord], image: str) -> list[float]:
    """Run the installed canonical image and load its scores.

    Returns:
        Computed scores in input order.

    Raises:
        RuntimeError: If the scoring process fails.
    """
    with tempfile.TemporaryDirectory(prefix="azimuth-docker-") as temp_dir:
        temp_root = Path(temp_dir)
        input_path = temp_root / "input.json"
        output_path = temp_root / "output.json"
        runner_path = temp_root / "runner.py"
        write_record_batch(input_path, records)
        write_helper_script(runner_path, DOCKER_RUNNER)
        command = build_docker_command(
            image=image,
            work_dir=temp_root,
        )
        result = subprocess.run(  # noqa: S603 - operator-selected programs, literal argv
            command, check=False, capture_output=True, text=True
        )
        if result.returncode != 0:
            message = (
                result.stderr.strip()
                or result.stdout.strip()
                or "docker backend failed"
            )
            raise RuntimeError(message)
        return json.loads(output_path.read_text(encoding="utf-8"))


def ensure_fork_checkout(cache_root: Path, url: str, ref: str) -> Path:
    """Reuse or atomically create the pinned source checkout.

    Returns:
        The pinned checkout directory.

    Raises:
        RuntimeError: If Git is absent or the revision differs.
    """
    git = shutil.which("git")
    if git is None:
        raise RuntimeError("git is required for the Azimuth fork fallback")
    checkout_dir = cache_root / ref
    if checkout_dir.exists():
        return checkout_dir
    cache_root.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix="fork-", dir=cache_root))
    try:
        subprocess.run(  # noqa: S603 - operator-selected programs, literal argv
            [git, "clone", url, str(temp_dir)],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(  # noqa: S603 - operator-selected programs, literal argv
            [git, "checkout", "--detach", ref],
            cwd=temp_dir,
            check=True,
            capture_output=True,
            text=True,
        )
        head = subprocess.run(  # noqa: S603 - operator-selected programs, literal argv
            [git, "rev-parse", "HEAD"],
            cwd=temp_dir,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        if head != ref:
            raise RuntimeError(f"Expected fork checkout {ref}, got {head}")
        temp_dir.replace(checkout_dir)
    except Exception:
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise
    return checkout_dir


def build_fork_command(
    *,
    helper_path: Path,
    checkout_dir: Path,
    input_path: Path,
    output_path: Path,
) -> list[str]:
    """Select the existing isolated Python 3.10 fork environment.

    Returns:
        Literal uv argv.
    """
    command = ["uv", "run", "--python", FORK_PYTHON]
    for dependency in FORK_RUNTIME_DEPENDENCIES:
        command.extend(["--with", dependency])
    command.extend(
        [
            "python",
            str(helper_path),
            str(checkout_dir),
            str(input_path),
            str(output_path),
        ]
    )
    return command


def run_fork_backend(
    records: list[SequenceRecord],
    *,
    cache_root: Path,
    url: str,
    ref: str,
) -> list[float]:
    """Run the pinned upstream fork and load its scores.

    Returns:
        Computed scores in input order.

    Raises:
        RuntimeError: If uv is absent or the scoring process fails.
    """
    if shutil.which("uv") is None:
        raise RuntimeError("uv is required for the Azimuth fork fallback")
    checkout_dir = ensure_fork_checkout(cache_root, url, ref)
    with tempfile.TemporaryDirectory(prefix="azimuth-fork-") as temp_dir:
        temp_root = Path(temp_dir)
        input_path = temp_root / "input.json"
        output_path = temp_root / "output.json"
        helper_path = temp_root / "run_fork.py"
        write_record_batch(input_path, records)
        write_helper_script(helper_path, FORK_RUNNER)
        command = build_fork_command(
            helper_path=helper_path,
            checkout_dir=checkout_dir,
            input_path=input_path,
            output_path=output_path,
        )
        env = dict(os.environ)
        env.setdefault("PYTHONWARNINGS", "ignore")
        result = subprocess.run(  # noqa: S603 - operator-selected programs, literal argv
            command, check=False, capture_output=True, text=True, env=env
        )
        if result.returncode != 0:
            message = (
                result.stderr.strip() or result.stdout.strip() or "fork backend failed"
            )
            raise RuntimeError(message)
        return json.loads(output_path.read_text(encoding="utf-8"))


def score_records(
    records: list[SequenceRecord],
    *,
    backend: Literal["docker", "fork"],
    image: str,
    cache_root: Path,
    fork_url: str,
    fork_ref: str,
) -> list[ScoredRecord]:
    """Attach upstream predictions to their original input records.

    Returns:
        Scores paired with their inputs and backend identity.
    """
    if backend == "docker":
        scores = run_docker_backend(records, image)
    else:
        scores = run_fork_backend(
            records, cache_root=cache_root, url=fork_url, ref=fork_ref
        )
    return [
        ScoredRecord(
            **record.model_dump(mode="json"),
            azimuth_score=float(score),
            backend=backend,
        )
        for record, score in zip(records, scores, strict=True)
    ]


def render_records(
    records: list[ScoredRecord], output_format: Literal["json", "csv", "tsv"]
) -> str:
    """Serialize scored records in the requested output format.

    Returns:
        Formatted result text.
    """
    rows = [record.model_dump(mode="json") for record in records]
    if output_format == "json":
        return json.dumps(rows, indent=2)
    delimiter = "," if output_format == "csv" else "\t"
    headers = [
        "record_id",
        "sequence",
        "aa_cut",
        "percent_peptide",
        "azimuth_score",
        "backend",
    ]
    lines = [delimiter.join(headers)]
    for row in rows:
        fields = [
            "" if row.get(header) is None else str(row[header]) for header in headers
        ]
        lines.append(delimiter.join(fields))
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    """Declare the established command-line arguments.

    Returns:
        The argument parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sequence",
        action="append",
        default=[],
        help="One Azimuth 30mer. Repeat for batches.",
    )
    parser.add_argument(
        "--aa-cut",
        type=float,
        default=None,
        help="Optional amino-acid cut position for one record.",
    )
    parser.add_argument(
        "--percent-peptide",
        type=float,
        default=None,
        help="Optional percent-peptide value for one record.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Optional JSON, JSONL, CSV, or TSV batch file.",
    )
    parser.add_argument(
        "--output-format", choices=("json", "csv", "tsv"), default="json"
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Optional output file. Defaults to stdout.",
    )
    parser.add_argument("--backend", choices=("auto", "docker", "fork"), default="auto")
    parser.add_argument("--docker-image", default=DEFAULT_DOCKER_IMAGE)
    parser.add_argument("--fork-url", default=DEFAULT_FORK_URL)
    parser.add_argument("--fork-ref", default=DEFAULT_FORK_REF)
    parser.add_argument(
        "--cache-root",
        type=Path,
        default=DEFAULT_CACHE_ROOT,
    )
    return parser


def main() -> None:
    """Run the command-line workflow with the caller-selected inputs."""
    parser = build_parser()
    args = parser.parse_args()
    try:
        records = load_records_from_args(args)
        backend = choose_backend(args.backend, args.docker_image)
        scored = score_records(
            records,
            backend=backend,
            image=args.docker_image,
            cache_root=args.cache_root,
            fork_url=args.fork_url,
            fork_ref=args.fork_ref,
        )
        payload = render_records(scored, args.output_format)
    except (
        RuntimeError,
        ValidationError,
        ValueError,
        subprocess.SubprocessError,
    ) as exc:
        parser.exit(1, f"{exc}\n")

    if args.out is None:
        sys.stdout.write(payload + "\n")
    else:
        args.out.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
