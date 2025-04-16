from io import StringIO
from pathlib import Path
from tempfile import NamedTemporaryFile

import pytest

from biov import BioDataFrame, read_gff3
from biov.io.gff import GFF_COLUMNS

# Test data
SAMPLE_GFF = """##gff-version 3
chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=TestGene;Dbxref=NCBI:123
chr1\t.\texon\t150\t180\t.\t+\t.\tID=exon1;Parent=gene1
"""


@pytest.fixture
def sample_gff():
    return read_gff3(StringIO(SAMPLE_GFF))


def test_read_gff3():
    df = read_gff3(StringIO(SAMPLE_GFF), explode_attributes=False)
    assert isinstance(df, BioDataFrame)
    assert list(df.columns) == GFF_COLUMNS
    assert len(df) == 2
    df = read_gff3(StringIO(SAMPLE_GFF), explode_attributes=True)
    assert "ID" in df.columns
    assert "Name" in df.columns
    assert "Dbxref" in df.columns
    assert "Parent" in df.columns


def test_to_gff3(sample_gff: BioDataFrame):
    with NamedTemporaryFile(mode="w", suffix=".gff") as f:
        sample_gff.to_gff3(f.name)
        content = Path(f.name).read_text()
        assert content.startswith("##gff-version 3")
        assert "chr1\t.\tgene\t100\t200" in content
    assert sample_gff.to_gff3() == SAMPLE_GFF
    sample_gff["hello"] = "world"
    assert (
        sample_gff.to_gff3()
        == """##gff-version 3
chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=TestGene;Dbxref=NCBI:123;hello=world
chr1\t.\texon\t150\t180\t.\t+\t.\tID=exon1;Parent=gene1;hello=world
"""
    )
