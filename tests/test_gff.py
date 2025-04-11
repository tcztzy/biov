from io import StringIO
from pathlib import Path
from tempfile import NamedTemporaryFile

import pytest

from biov.gff import GFF3_COLUMNS, GFFDataFrame, read_gff3

# Test data
SAMPLE_GFF = """##gff-version 3
chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=TestGene;Dbxref=NCBI:123
chr1\t.\texon\t150\t180\t.\t+\t.\tID=exon1;Parent=gene1
"""


@pytest.fixture
def sample_gff_df():
    data = {
        "seqid": ["chr1", "chr1"],
        "source": [".", "."],
        "type": ["gene", "exon"],
        "start": [100, 150],
        "end": [200, 180],
        "score": [".", "."],
        "strand": ["+", "+"],
        "phase": [".", "."],
        "attributes": [
            "ID=gene1;Name=TestGene;Dbxref=NCBI:123",
            "ID=exon1;Parent=gene1",
        ],
    }
    return GFFDataFrame(data)


def test_read_gff3():
    df = read_gff3(StringIO(SAMPLE_GFF))
    assert isinstance(df, GFFDataFrame)
    assert list(df.columns) == GFF3_COLUMNS
    assert len(df) == 2


def test_gffdataframe_validation():
    # Should pass with all required columns
    valid_data = {col: [] for col in GFF3_COLUMNS}
    GFFDataFrame(valid_data)

    # Should fail if missing columns
    invalid_data = {col: [] for col in GFF3_COLUMNS[:-1]}
    with pytest.raises(AttributeError):
        GFFDataFrame(invalid_data)


def test_to_gff3(sample_gff_df):
    with NamedTemporaryFile(mode="w", suffix=".gff") as f:
        sample_gff_df.to_gff3(f.name)
        content = Path(f.name).read_text()
        assert content.startswith("##gff-version 3\n")
        assert "chr1\t.\tgene\t100\t200" in content


def test_attributes_to_columns(sample_gff_df):
    df = sample_gff_df.attributes
    assert "ID" in df.columns
    assert "Name" in df.columns
    assert "Parent" in df.columns
    assert df.loc[0, "ID"] == "gene1"
    assert df.loc[0, "Name"] == "TestGene"
    assert df.loc[1, "Parent"] == "gene1"
