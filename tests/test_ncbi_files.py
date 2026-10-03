"""NCBI file downloads retain full data and fail on unavailable representations."""

import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from test_artifacts import Response

import biov.ncbi_files as ncbi


@pytest.mark.parametrize(
    ("namespace", "accession", "kind", "payload", "url_fragment"),
    [
        (
            "clinvar",
            "9",
            "variant_xml",
            (
                b'<ClinVarResult-Set><VariationArchive VariationID="9">'
                b"<ClassifiedRecord/></VariationArchive></ClinVarResult-Set>"
            ),
            "is_variationid=true",
        ),
        (
            "dbsnp",
            "rs328",
            "refsnp_json",
            b'{"refsnp_id":"328","primary_snapshot_data":{"placements_with_allele":[]}}',
            "/variation/v0/refsnp/328",
        ),
        (
            "geo",
            "GSM575",
            "soft",
            b"^SAMPLE = GSM575\n!sample_table_begin\n",
            "targ=self",
        ),
        (
            "geo",
            "GPL96",
            "soft",
            b"^PLATFORM = GPL96\n!platform_table_begin\n",
            "targ=self",
        ),
        (
            "geo",
            "GDS100",
            "soft",
            b"^DATASET = GDS100\n",
            "GDSnnn/GDS100/soft/GDS100_full.soft.gz",
        ),
        (
            "geo",
            "GSE100",
            "soft",
            b"^SERIES = GSE100\n",
            "GSEnnn/GSE100/soft/GSE100_family.soft.gz",
        ),
    ],
)
def test_download_complete_ncbi_records(
    monkeypatch, tmp_path, namespace, accession, kind, payload, url_fragment
):
    """Keep full responses and decompress SOFT without replacing their content."""
    calls = []

    def open_url(url, *, timeout):
        calls.append(url)
        return Response(gzip.compress(payload) if url.endswith(".gz") else payload)

    monkeypatch.setattr(ncbi, "urlopen", open_url)
    target = tmp_path / "data"
    ncbi.download_ncbi_file(namespace, accession, kind, target)
    assert target.read_bytes() == payload
    assert len(calls) == 1 and url_fragment in calls[0]


def test_pubmed_uses_current_cloud_version_and_verifies_pdf(monkeypatch, tmp_path):
    """Use the current version, exact PMID, PDF link, and published checksum."""
    pdf = b"%PDF-1.4\narticle content\n%%EOF\n"
    checksum = hashlib.md5(pdf, usedforsecurity=False).hexdigest()
    metadata = {
        "pmcid": "PMC3531190",
        "version": 2,
        "pmid": 23193287,
        "pdf_url": "s3://pmc-oa-opendata/PMC3531190.2/PMC3531190.2.pdf?md5=" + checksum,
    }
    calls = []

    def open_url(url, *, timeout):
        calls.append(url)
        if "idconv" in url:
            result = {
                "records": [
                    {
                        "pmid": 23193287,
                        "versions": [
                            {"pmcid": "PMC3531190.1", "current": False},
                            {"pmcid": "PMC3531190.2", "current": True},
                        ],
                    }
                ]
            }
        elif "/metadata/" in url:
            assert url.endswith("PMC3531190.2.json")
            result = metadata
        else:
            assert url.endswith("PMC3531190.2/PMC3531190.2.pdf")
            return Response(pdf)
        return Response(json.dumps(result).encode())

    monkeypatch.setattr(ncbi, "urlopen", open_url)
    target = tmp_path / "article.pdf"
    ncbi.download_ncbi_file("pubmed", "23193287", "article_pdf", target)
    assert target.read_bytes() == pdf
    assert len(calls) == 3
    metadata["pdf_url"] = None
    with pytest.raises(FileNotFoundError, match="no downloadable PMC PDF"):
        ncbi.download_ncbi_file(
            "pubmed", "23193287", "article_pdf", tmp_path / "absent"
        )


def test_geo_requires_one_nonempty_matching_matrix(monkeypatch, tmp_path):
    """Keep table metadata and reject multiple, missing, empty, or mismatched data."""
    directory = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE100/matrix/"
    matrix_url = directory + "GSE100_series_matrix.txt.gz"
    candidates = [matrix_url]
    matrix = (
        b'!Series_geo_accession\t"GSE100"\n!series_matrix_table_begin\n'
        b'"ID_REF"\t"GSM575"\n"a"\t1.2\n!series_matrix_table_end\n'
    )
    monkeypatch.setattr(
        ncbi.fsspec,
        "filesystem",
        lambda _: SimpleNamespace(ls=lambda url, detail: candidates),
    )
    monkeypatch.setattr(
        ncbi, "urlopen", lambda url, timeout: Response(gzip.compress(matrix))
    )
    target = tmp_path / "matrix.txt"
    ncbi.download_ncbi_file("geo", "GSE100", "expression_matrix", target)
    assert target.read_bytes() == matrix
    for urls in ([], [matrix_url, directory + "GSE100-GPL96_series_matrix.txt.gz"]):
        candidates[:] = urls
        with pytest.raises(FileNotFoundError, match="Series matrices"):
            ncbi.download_ncbi_file(
                "geo", "GSE100", "expression_matrix", tmp_path / "other"
            )
    with pytest.raises(FileNotFoundError, match="artifact='soft'"):
        ncbi.download_ncbi_file(
            "geo", "GSM575", "expression_matrix", tmp_path / "sample"
        )
    for invalid in (
        matrix.replace(b'"GSE100"', b'"GSE200"'),
        matrix.replace(b'"a"\t1.2\n', b""),
    ):
        target.write_bytes(invalid)
        with pytest.raises(ValueError, match="matching expression_matrix"):
            ncbi.validate_ncbi_file("geo", "GSE100", "expression_matrix", target)


@pytest.mark.parametrize(
    ("namespace", "accession", "kind", "payload"),
    [
        ("pubmed", "1", "article_pdf", b"<html>Error</html>"),
        (
            "clinvar",
            "9",
            "variant_xml",
            b'<ClinVarResult-Set><VariationArchive VariationID="10"/></ClinVarResult-Set>',
        ),
        ("dbsnp", "rs328", "refsnp_json", b'{"refsnp_id":"1"}'),
        ("geo", "GSM575", "soft", b"<html>Error</html>"),
    ],
)
def test_ncbi_cache_validation_rejects_wrong_record(
    tmp_path: Path, namespace, accession, kind, payload
):
    """Prevent an HTML error or another accession becoming a cached artifact."""
    target = tmp_path / "bad"
    target.write_bytes(payload)
    with pytest.raises(ValueError):
        ncbi.validate_ncbi_file(namespace, accession, kind, target)
