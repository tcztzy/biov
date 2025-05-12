"""Tests for range operations in BioDataFrame."""

# ruff: noqa: S101,DOC201,PD901

import pyranges as pr
import pytest

from biov import BioDataFrame
from biov.ranges import from_pyranges, to_pyranges


@pytest.fixture
def sample_ranges():
    """Create a sample BioDataFrame with genomic ranges."""
    data = {
        "seqid": ["chr1", "chr1", "chr2", "chr2", "chr3"],
        "start": [100, 200, 150, 500, 300],
        "end": [150, 250, 200, 550, 350],
        "strand": ["+", "-", "+", "+", "-"],
        "name": ["A", "B", "C", "D", "E"],
        "score": [10, 20, 15, 25, 30],
    }
    return BioDataFrame(data)


@pytest.fixture
def other_ranges():
    """Create another sample BioDataFrame with genomic ranges."""
    data = {
        "seqid": ["chr1", "chr1", "chr2", "chr3", "chr3"],
        "start": [120, 250, 180, 320, 400],
        "end": [170, 300, 220, 380, 450],
        "strand": ["+", "+", "-", "-", "+"],
        "name": ["X", "Y", "Z", "W", "V"],
        "score": [5, 15, 10, 20, 25],
    }
    return BioDataFrame(data)


def test_to_pyranges(sample_ranges):
    """Test conversion from BioDataFrame to PyRanges."""
    gr = to_pyranges(sample_ranges)
    assert isinstance(gr, pr.PyRanges)
    assert len(gr) == len(sample_ranges)
    assert "Chromosome" in gr.as_df().columns
    assert "Start" in gr.as_df().columns
    assert "End" in gr.as_df().columns
    assert "Strand" in gr.as_df().columns
    assert "name" in gr.as_df().columns
    assert "score" in gr.as_df().columns


def test_from_pyranges(sample_ranges):
    """Test conversion from PyRanges to BioDataFrame."""
    gr = to_pyranges(sample_ranges)
    df = from_pyranges(gr)
    assert isinstance(df, BioDataFrame)
    assert len(df) == len(sample_ranges)
    assert "seqid" in df.columns
    assert "start" in df.columns
    assert "end" in df.columns
    assert "strand" in df.columns
    assert "name" in df.columns
    assert "score" in df.columns


def test_biodataframe_to_pyranges(sample_ranges):
    """Test BioDataFrame.to_pyranges method."""
    gr = sample_ranges.to_pyranges()
    assert isinstance(gr, pr.PyRanges)
    assert len(gr) == len(sample_ranges)


def test_overlap(sample_ranges, other_ranges):
    """Test overlap operation."""
    result = sample_ranges.overlap(other_ranges)
    assert isinstance(result, BioDataFrame)
    # There should be at least one overlap between the two sets
    assert len(result) > 0
    # Check that the result has columns from both DataFrames
    assert "name" in result.columns


def test_intersect(sample_ranges, other_ranges):
    """Test intersect operation."""
    result = sample_ranges.intersect(other_ranges)
    assert isinstance(result, BioDataFrame)
    # The result should contain only the overlapping regions
    for _, row in result.iterrows():
        # For each result row, check that it's a valid intersection
        # by ensuring it exists within the original ranges
        assert any(
            (
                row["seqid"] == s_row["seqid"]
                and row["start"] >= s_row["start"]
                and row["end"] <= s_row["end"]
            )
            for _, s_row in sample_ranges.iterrows()
        )


def test_subtract_ranges(sample_ranges, other_ranges):
    """Test subtract_ranges operation."""
    result = sample_ranges.subtract_ranges(other_ranges)
    assert isinstance(result, BioDataFrame)
    # The result should not contain any regions that fully overlap with other_ranges
    for _, row in result.iterrows():
        # For each result row, check that it doesn't fully overlap with any row in other_ranges
        assert not any(
            (
                row["seqid"] == o_row["seqid"]
                and row["start"] >= o_row["start"]
                and row["end"] <= o_row["end"]
            )
            for _, o_row in other_ranges.iterrows()
        )


def test_nearest(sample_ranges, other_ranges):
    """Test nearest operation."""
    result = sample_ranges.nearest(other_ranges)
    assert isinstance(result, BioDataFrame)
    # PyRanges nearest may not return a result for every input range
    # if there's no nearby range, so we don't check the length
    # Check that the result has columns from both DataFrames
    assert "name" in result.columns
    assert "name_b" in result.columns


def test_custom_column_names():
    """Test using custom column names."""
    data = {
        "chrom": ["chr1", "chr2"],
        "chromStart": [100, 200],
        "chromEnd": [150, 250],
        "name": ["A", "B"],
    }
    df = BioDataFrame(data)

    # Test conversion with custom column names
    gr = to_pyranges(df, seqid_col="chrom", start_col="chromStart", end_col="chromEnd")
    assert isinstance(gr, pr.PyRanges)
    assert len(gr) == len(df)

    # Test conversion back
    result = from_pyranges(
        gr, seqid_col="chrom", start_col="chromStart", end_col="chromEnd"
    )
    assert isinstance(result, BioDataFrame)
    assert "chrom" in result.columns
    assert "chromStart" in result.columns
    assert "chromEnd" in result.columns
