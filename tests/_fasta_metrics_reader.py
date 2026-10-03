"""Read moved native FASTA metric bundles with stock PyArrow and no BioV."""

import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
import re
import sys
from fractions import Fraction
from pathlib import Path

assert os.environ["PATH"] == ""
assert os.environ["PYTHONPATH"] == ""
assert sys.prefix != sys.base_prefix, "Reader must use an independent virtualenv"
assert importlib.util.find_spec("biov") is None, "BioV must not be installed"
assert {
    distribution.metadata["Name"].lower()
    for distribution in importlib.metadata.distributions()
} == {"pyarrow"}, "Provision a clean pyarrow-only virtualenv"


def offline_only(event, args):
    """Reject network, process, database and BioV dependencies in the reader.

    Raises:
        AssertionError: If the independent reader attempts a forbidden operation.
    """
    if event.startswith(("socket.", "subprocess.", "sqlite3.")) or event in {
        "os.system",
        "os.exec",
        "os.posix_spawn",
        "os.fork",
        "os.forkpty",
    }:
        raise AssertionError(f"Standalone reading must be offline: {event}")
    if event == "import" and args[0].split(".")[0] in {"biov", "sqlite3"}:
        raise AssertionError(f"Forbidden reader dependency: {args[0]}")


sys.addaudithook(offline_only)
import pyarrow as pa
import pyarrow.compute as pc
from pyarrow import ipc

assert pa.__version__ == "25.0.1", pa.__version__
COMPUTE = {
    name: getattr(pc, name) for name in ("equal", "sort_indices", "sum", "is_null")
}
COLUMNS = [
    "sequence_id",
    "start",
    "end",
    "length",
    "is_full_window",
    "canonical_base_count",
    "gc_base_count",
    "gc_fraction",
    "weighted_gc_fraction",
]
DTYPES = ["str", "i64", "i64", "i64", "bool", "i64", "i64", "f64", "f64"]
ARROW_TYPES = [
    pa.large_string(),
    pa.int64(),
    pa.int64(),
    pa.int64(),
    pa.bool_(),
    pa.int64(),
    pa.int64(),
    pa.float64(),
    pa.float64(),
]
IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def expected_synthetic():
    """Return independent complete fixture records without a producer preview."""
    sequence_id = "chrAlpha.1"
    sequence = "ACGTacgtnNRYacGTACGTnnnnACGTacGTACGTA"
    expected = []
    for start in range(0, len(sequence), 4):
        bases = sequence[start : start + 4].upper()
        canonical = sum(base in "ACGT" for base in bases)
        gc = sum(base in "GC" for base in bases)
        weighted = sum(
            (
                Fraction(sum(base in "GC" for base in IUPAC[code]), len(IUPAC[code]))
                for code in bases
            ),
            start=Fraction(),
        )
        expected.append(
            {
                "sequence_id": sequence_id,
                "start": start,
                "end": start + len(bases),
                "length": len(bases),
                "is_full_window": len(bases) == 4,
                "canonical_base_count": canonical,
                "gc_base_count": gc,
                "gc_fraction": gc / canonical if canonical else None,
                "weighted_gc_fraction": float(weighted / len(bases)),
            }
        )
    return expected


def assert_rows(actual, expected):
    """Check complete typed rows with a narrow float tolerance."""
    assert len(actual) == len(expected)
    for observed, reference in zip(actual, expected, strict=True):
        assert observed.keys() == reference.keys()
        for key, value in reference.items():
            if isinstance(value, float):
                assert math.isclose(observed[key], value, rel_tol=1e-14, abs_tol=1e-14)
            else:
                assert observed[key] == value, (key, observed, reference)


def read_bundle(root):
    """Verify, discover and completely read one four-file metric artifact.

    Returns:
        The complete Arrow table, strict record and descriptive manifest.
    """
    manifests = list(root.glob("artifact_*.manifest.json"))
    assert len(manifests) == 1
    manifest_path = manifests[0]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["manifest_version"] == 1
    assert manifest["format"] == "arrow_ipc_file"
    assert set(manifest["files"]) == {"arrow", "record", "readme"}
    paths = {}
    for role, filename in manifest["files"].items():
        relative = Path(filename)
        assert relative.parts == (filename,) and "\\" not in filename
        assert not relative.is_absolute() and filename not in {".", ".."}
        path = root / filename
        assert path.is_file() and not path.is_symlink()
        paths[role] = path
    assert set(root.iterdir()) == {manifest_path, *paths.values()}
    arrow_bytes, record_bytes = (
        paths["arrow"].read_bytes(),
        paths["record"].read_bytes(),
    )
    assert len(arrow_bytes) == manifest["content"]["bytes"]
    assert hashlib.sha256(arrow_bytes).hexdigest() == manifest["content"]["sha256"]
    assert len(record_bytes) == manifest["record"]["bytes"]
    assert hashlib.sha256(record_bytes).hexdigest() == manifest["record"]["sha256"]
    record = json.loads(record_bytes)
    assert record["record_version"] == 3
    assert record["sha256"] == manifest["content"]["sha256"]
    assert record["bytes"] == len(arrow_bytes)
    assert record["row_count"] == manifest["content"]["row_count"]
    assert record["schema"] == [
        {"name": name, "dtype": dtype}
        for name, dtype in zip(COLUMNS, DTYPES, strict=True)
    ]
    assert "csv_schema_policy" not in record["provenance"]
    assert record["provenance"]["declared_schema"] == {}
    origin = record["provenance"]["sequence_origin"]
    assert origin == manifest["lineage"]["sequence_origin"]
    assert origin["format"] == manifest["source"]["format"] == "fasta"
    assert origin["reference"] == "refseq.gcf:GCF_000005845.2"
    assert origin["snapshot_id"].startswith("sha256-")
    assert origin["recipe_id"].startswith("sha256-")
    assert origin["recipe_id"] != origin["snapshot_id"]
    for field in ("fai_sha256", "dictionary_sha256"):
        assert re.fullmatch("[a-f0-9]{64}", origin[field]), (field, origin)
    assert origin["coordinates"] == "0-based-half-open;source-sequence-relative"
    assert origin["units"] == "length/counts:bases;gc-fractions:dimensionless"
    assert origin["canonical_gc_policy"] == (
        "(G+C)/(A+C+G+T);ascii-case-insensitive;ambiguity-excluded;zero-denominator-null"
    )
    assert origin["weighted_gc_policy"] == (
        "biov-core-iupac-dna-gc;all-bases-denominator;N=1/2;B,V=2/3;D,H=1/3"
    )
    assert origin["algorithm"] == "nonoverlapping-fasta-gc-windows"
    assert origin["algorithm_revision"] == 1
    assert manifest["source"]["required_for_reading"] is False
    assert manifest["source"]["sha256"] == record["provenance"]["source_sha256"]
    assert not (root / manifest["source"]["historical_path"]).exists()
    assert manifest["lineage"]["operations"] == record["provenance"]["operations"]
    assert manifest["lineage"]["historical_references_only"] is True
    assert manifest["scientific_metadata"] == record["scientific_metadata"]
    assert record["scientific_metadata"]["species"] is None
    assert manifest["versions"] == {"reference_version": None, "provider_release": None}
    readme = paths["readme"].read_text(encoding="utf-8")
    examples = re.findall(
        r"^```python\s*\n(.*?)^```\s*$", readme, re.MULTILINE | re.DOTALL
    )
    assert len(examples) == 1
    previous = Path.cwd()
    os.chdir(root)
    try:
        exec(  # noqa: S102
            compile(examples[0], f"{paths['readme']} [Python quick-start]", "exec"),
            {"__name__": "__main__"},
        )
    finally:
        os.chdir(previous)
    with pa.BufferReader(arrow_bytes) as reader:
        table = ipc.open_file(reader).read_all()
    assert table.column_names == COLUMNS
    assert table.num_rows == record["row_count"]
    assert len(manifest["columns"]) == len(COLUMNS)
    for field, expected, column in zip(
        table.schema, ARROW_TYPES, manifest["columns"], strict=True
    ):
        assert field.type == expected, (field, expected)
        assert field.name == column["name"] and field.nullable
        assert table[field.name].null_count == column["null_count"]
        assert column["description"], (
            "Metric columns need their known scientific meanings"
        )
        if field.name in {"start", "end"}:
            assert column["coordinates"] == origin["coordinates"]
        if field.name in {
            "start",
            "end",
            "length",
            "canonical_base_count",
            "gc_base_count",
        }:
            assert column["units"] == "bases"
        if field.name in {"gc_fraction", "weighted_gc_fraction"}:
            assert column["units"] == "dimensionless"
    return table, record, manifest


mode = sys.argv[1]
assert mode in {"synthetic", "real"}
root = Path.cwd()
full, full_record, full_manifest = read_bundle(root / "full")
selected, selected_record, _ = read_bundle(root / "selected")
reopened, reopened_record, _ = read_bundle(root / "reopened")
assert (
    selected_record["provenance"]["sequence_origin"]
    == full_record["provenance"]["sequence_origin"]
)
assert (
    reopened_record["provenance"]["sequence_origin"]
    == full_record["provenance"]["sequence_origin"]
)
assert (
    reopened_record["reopen_verification"]["original_provenance"]
    == "recorded_claims_not_independently_verified"
)
assert selected_record["provenance"]["operations"]
assert reopened_record["provenance"]["operations"]
assert_rows(reopened.to_pylist(), selected.to_pylist())
# All comparisons use stock Arrow kernels directly, without compatibility casts.
kept = full.filter(COMPUTE["equal"](full["is_full_window"], pa.scalar(True)))
indices = COMPUTE["sort_indices"](kept, sort_keys=[("start", "descending")])
assert_rows(kept.take(indices).to_pylist(), selected.to_pylist())
assert (
    COMPUTE["sum"](full["length"]).as_py()
    == full_record["provenance"]["sequence_origin"]["sequence_length"]
)
assert COMPUTE["sum"](selected["is_full_window"]).as_py() == selected.num_rows
if mode == "synthetic":
    expected = expected_synthetic()
    assert_rows(full.to_pylist(), expected)
    assert_rows(
        selected.to_pylist(), [row for row in expected if row["is_full_window"]][::-1]
    )
    assert full.num_rows > 5, "Analyze records beyond a default preview"
    ambiguous = full.filter(COMPUTE["is_null"](full["gc_fraction"]))
    assert ambiguous.num_rows == 2
    assert all(
        row["canonical_base_count"] == row["gc_base_count"] == 0
        and math.isclose(row["weighted_gc_fraction"], 0.5, rel_tol=1e-14)
        for row in ambiguous.to_pylist()
    )
    assert COMPUTE["sum"](full["canonical_base_count"]).as_py() == 29
    assert COMPUTE["sum"](full["gc_base_count"]).as_py() == 14
else:
    assert full.num_rows == 465 and selected.num_rows == 464
    for index, row in enumerate(full.to_pylist()):
        assert row["sequence_id"] == "NC_000913.3"
        assert row["start"] == index * 10_000
        assert row["end"] == min((index + 1) * 10_000, 4_641_652)
        assert row["length"] == row["end"] - row["start"]
        assert row["is_full_window"] is (index < 464)
        assert row["canonical_base_count"] == row["length"]
        assert math.isclose(
            row["gc_fraction"], row["gc_base_count"] / row["length"], rel_tol=1e-14
        )
    assert COMPUTE["sum"](full["length"]).as_py() == 4_641_652
    assert COMPUTE["sum"](full["gc_base_count"]).as_py() == 2_357_528
    assert COMPUTE["sum"](full["canonical_base_count"]).as_py() == 4_641_652
    assert full["start"][-1].as_py() == 4_640_000
    assert full["end"][-1].as_py() == 4_641_652
    assert full["length"][-1].as_py() == 1_652
    assert full["is_full_window"][-1].as_py() is False
    assert full["gc_fraction"].null_count == 0
    assert all(
        row["gc_fraction"] == row["weighted_gc_fraction"] for row in full.to_pylist()
    )
assert not any(name == "biov" or name.startswith("biov.") for name in sys.modules)
assert "sqlite3" not in sys.modules
print("standalone FASTA metrics filter and summary passed")  # noqa: T201
