"""Acceptance tests for BioV genomic range semantics."""

import importlib
import tomllib
from collections.abc import Sequence
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import biov
from biov import BioDataFrame
from biov.ranges import NearestHow


def _ranges(**columns: Sequence[object] | pd.Series) -> BioDataFrame:
    return BioDataFrame(columns)


def test_overlap_uses_half_open_coordinates_strand_groups_and_stable_order() -> None:
    """Select each matching input row once in its original order."""
    query = _ranges(
        seqid=["chr1", "chr1", "chr2", "chr1", "chr1"],
        start=[0, 0, 0, 0, 20],
        end=[10, 10, 10, 10, 30],
        strand=["+", "-", "+", "+", "+"],
        name=["plus", "minus", "other-chromosome", "duplicate", "touch-only"],
    )
    other = _ranges(
        seqid=["chr1", "chr1", "chr2", "chr1"],
        start=[10, 5, 2, 4],
        end=[15, 6, 3, 5],
        strand=["+", "-", "+", "+"],
    )

    result = query.overlap(other)

    assert result["name"].tolist() == [
        "plus",
        "minus",
        "other-chromosome",
        "duplicate",
    ]
    assert isinstance(result, BioDataFrame)
    assert isinstance(result.index, pd.RangeIndex)


def test_overlap_modes_have_explicit_containment_direction() -> None:
    """Distinguish self-containing from self-contained overlap modes."""
    query = _ranges(
        seqid=["chr1", "chr1", "chr1"],
        start=[0, 5, 40],
        end=[20, 10, 50],
        name=["contains", "member", "none"],
    )
    other = _ranges(seqid=["chr1", "chr1"], start=[2, 0], end=[15, 30])

    assert query.overlap(other, how="first")["name"].tolist() == [
        "contains",
        "member",
    ]
    assert query.overlap(other, how="last")["name"].tolist() == [
        "contains",
        "member",
    ]
    assert query.overlap(other, how="containment")["name"].tolist() == ["contains"]
    assert query.overlap(other, how="member")["name"].tolist() == [
        "contains",
        "member",
    ]


def test_intersect_expands_pairs_in_input_order_and_preserves_duplicates() -> None:
    """Emit one clipped self row for every overlap pair."""
    query = _ranges(
        seqid=["chr1", "chr1", "chr1"],
        start=[10, 0, 0],
        end=[20, 10, 10],
        name=["late", "early", "early-duplicate"],
    )
    other = _ranges(seqid=["chr1", "chr1"], start=[8, 2], end=[12, 4])

    result = query.intersect(other)

    assert result[["name", "start", "end"]].to_dict(  # pyright: ignore[reportCallIssue]
        orient="records"
    ) == [
        {"name": "late", "start": 10, "end": 12},
        {"name": "early", "start": 8, "end": 10},
        {"name": "early", "start": 2, "end": 4},
        {"name": "early-duplicate", "start": 8, "end": 10},
        {"name": "early-duplicate", "start": 2, "end": 4},
    ]


def test_subtract_ranges_unions_masks_and_orders_duplicate_fragments() -> None:
    """Subtract unioned masks without collapsing duplicate query rows."""
    query = _ranges(
        seqid=["chr1", "chr1", "chr1"],
        start=[10, 0, 0],
        end=[20, 10, 10],
        name=["late", "early", "early-duplicate"],
    )
    other = _ranges(
        seqid=["chr1", "chr1", "chr1", "chr1"],
        start=[5, 3, 5, 12],
        end=[7, 6, 7, 15],
    )

    result = query.subtract_ranges(other)

    assert result[["name", "start", "end"]].to_dict(  # pyright: ignore[reportCallIssue]
        orient="records"
    ) == [
        {"name": "late", "start": 10, "end": 12},
        {"name": "late", "start": 15, "end": 20},
        {"name": "early", "start": 0, "end": 3},
        {"name": "early", "start": 7, "end": 10},
        {"name": "early-duplicate", "start": 0, "end": 3},
        {"name": "early-duplicate", "start": 7, "end": 10},
    ]


@pytest.mark.parametrize(
    ("how", "expected"),
    [
        ("next", ["plus-right-first", "minus-right", "overlap"]),
        ("previous", ["plus-left", "minus-left", "overlap"]),
        ("upstream", ["plus-left", "minus-right", "overlap"]),
        ("downstream", ["plus-right-first", "minus-left", "overlap"]),
    ],
)
def test_nearest_direction_ties_overlap_distance_and_missing_groups(
    how: NearestHow, expected: list[str]
) -> None:
    """Map public directions explicitly and choose the lowest input-row tie."""
    query = _ranges(
        seqid=["chr1", "chr1", "chr1", "chr2"],
        start=[10, 10, 40, 0],
        end=[20, 20, 50, 10],
        strand=["+", "-", "+", "+"],
        name=["plus", "minus", "query-overlap", "no-candidate"],
    )
    other = _ranges(
        seqid=["chr1"] * 6,
        start=[20, 20, 5, 20, 5, 45],
        end=[25, 25, 10, 25, 10, 46],
        strand=["+", "+", "+", "-", "-", "+"],
        name=[
            "plus-right-first",
            "plus-right-second",
            "plus-left",
            "minus-right",
            "minus-left",
            "overlap",
        ],
        score=[0, 1, 2, 3, 4, 5],
    )

    result = query.nearest(other, how=how)

    assert result["name_b"].tolist() == expected
    assert result["Distance"].tolist() == [1, 1, 0]
    assert result["score"].tolist() == [
        other.loc[other["name"].eq(name), "score"].iloc[0] for name in expected
    ]


def test_custom_columns_and_explicit_unstranded_grouping() -> None:
    """Support custom coordinate names and explicit strand ignoring."""
    query = _ranges(chrom=["chr1"], chromStart=[0], chromEnd=[10], orientation=["+"])
    other = _ranges(chrom=["chr1"], chromStart=[2], chromEnd=[3], orientation=["-"])

    result = query.intersect(
        other,
        seqid_col="chrom",
        start_col="chromStart",
        end_col="chromEnd",
        strand_col=None,
    )

    assert result[["chromStart", "chromEnd"]].to_numpy().tolist() == [[2, 3]]


@pytest.mark.parametrize(
    ("query", "other", "message"),
    [
        (
            _ranges(seqid=["chr1"], start=[0], end=[1], strand=["+"]),
            _ranges(seqid=["chr1"], start=[0], end=[1]),
            "strand column",
        ),
        (
            _ranges(seqid=["chr1"], start=[0], end=[1], strand=["."]),
            _ranges(seqid=["chr1"], start=[0], end=[1], strand=["."]),
            "strand",
        ),
        (
            _ranges(seqid=["chr1"], start=[1], end=[1]),
            _ranges(seqid=["chr1"], start=[0], end=[1]),
            "start < end",
        ),
    ],
)
def test_range_validation_is_stable(
    query: BioDataFrame, other: BioDataFrame, message: str
) -> None:
    """Reject ambiguous grouping and malformed half-open intervals."""
    with pytest.raises(ValueError, match=message):
        query.overlap(other)


def test_empty_range_inputs_have_stable_schemas() -> None:
    """Handle empty inputs without entering an interval kernel."""
    query = _ranges(seqid=["chr1"], start=[0], end=[10], name=["q"])
    empty_query = query.iloc[:0]
    other = _ranges(seqid=["chr1"], start=[2], end=[3], name=["r"], score=[1])
    empty_other = other.iloc[:0]

    assert_frame_equal(query.overlap(empty_other), query.iloc[:0])
    assert_frame_equal(query.intersect(empty_other), query.iloc[:0])
    assert_frame_equal(query.subtract_ranges(empty_other), query.reset_index(drop=True))
    assert empty_query.overlap(other).empty
    nearest = query.nearest(empty_other)
    assert nearest.empty
    assert nearest.columns.tolist() == [
        "seqid",
        "start",
        "end",
        "name",
        "start_b",
        "end_b",
        "name_b",
        "score",
        "Distance",
    ]

    typed_query = _ranges(
        seqid=pd.Series(["chr1"], dtype="string"),
        start=pd.Series([0], dtype="Int32"),
        end=pd.Series([10], dtype="Int32"),
    )
    typed_empty_other = typed_query.iloc[:0]
    assert_frame_equal(
        typed_query.intersect(typed_empty_other),
        typed_query.iloc[:0],
    )


def test_pyranges_compatibility_surface_is_removed() -> None:
    """Expose no PyRanges conversion or object-accepting compatibility path."""
    ranges_module = importlib.import_module("biov.ranges")
    query = _ranges(seqid=["chr1"], start=[0], end=[1])

    assert not hasattr(biov, "to_pyranges")
    assert not hasattr(biov, "from_pyranges")
    assert not hasattr(ranges_module, "to_pyranges")
    assert not hasattr(ranges_module, "from_pyranges")
    assert not hasattr(query, "to_pyranges")
    with pytest.raises(TypeError, match="BioDataFrame"):
        query.overlap(object())  # type: ignore[arg-type]


def test_legacy_interval_dependencies_are_absent() -> None:
    """Keep PyRanges 0.x and its implementation dependencies out of the lock."""
    root = Path(__file__).parent.parent
    project = tomllib.loads((root / "pyproject.toml").read_text())
    lock = tomllib.loads((root / "uv.lock").read_text())
    forbidden = {"pyranges", "sorted-nearest", "ncls", "setuptools"}
    direct = {
        dependency.split("[", maxsplit=1)[0].split(">", maxsplit=1)[0]
        for dependency in project["project"]["dependencies"]
    }
    locked = {package["name"] for package in lock["package"]}

    assert project["project"]["requires-python"] == ">=3.12"
    assert forbidden.isdisjoint(direct)
    assert forbidden.isdisjoint(locked)
