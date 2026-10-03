"""Native analysis files, identifier validation, and shared cache behavior."""

import hashlib
import json

import pytest
from test_artifacts import Response

from biov import artifacts, file_sources, ncbi_files, path


@pytest.mark.parametrize(
    ("uri", "kind", "payload"),
    [
        (
            "pdb://1CRN",
            "structure_cif",
            b"data_1CRN\n_entry.id 1CRN\n_atom_site.Cartn_x 1.2\n",
        ),
        ("ensembl://ENSG00000139618", "sequence_fasta", b">ENSG00000139618.19\nACGT\n"),
        (
            "chembl.compound://CHEMBL25",
            "structure_sdf",
            b"\nM  END\n> <chembl_id>\nCHEMBL25\n\n$$$$\n",
        ),
        ("pubchem.compound://2244", "structure_sdf", b"2244\nM  END\n$$$$\n"),
        ("arxiv://hep-th/9901001", "article_pdf", b"%PDF-1.7\n"),
        ("emdb://EMD-1001", "density_map", b"\0" * 208 + b"MAP " + b"\0" * 812),
        (
            "clinicaltrials://NCT00222573",
            "study_json",
            json.dumps(
                {"protocolSection": {"identificationModule": {"nctId": "NCT00222573"}}}
            ).encode(),
        ),
        (
            "dailymed://973a9333-fec7-46dd-8eb5-25738f06ee54",
            "label_xml",
            b'<document xmlns="urn:hl7-org:v3"><setId root="973a9333-fec7-46dd-8eb5-25738f06ee54"/></document>',
        ),
        (
            "reactome://R-HSA-199420",
            "pathway_sbml",
            b'<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"><model name="R-HSA-199420"/></sbml>',
        ),
    ],
)
def test_native_files_validate_and_reuse_cache(
    uri, kind, payload, tmp_path, monkeypatch
):
    """Download, validate, and reuse one cached provider file per URI."""
    calls = []

    def download(namespace, accession, requested_kind, destination, properties):
        assert requested_kind == kind
        assert properties["url"].startswith("https://")
        destination.write_bytes(payload)
        calls.append(uri)

    monkeypatch.setattr(artifacts.settings, "home", tmp_path)
    monkeypatch.setattr(artifacts, "download_file", download)
    validations = []
    validate = artifacts.validate_file

    def checked(namespace, accession, selected_kind, source):
        validations.append(selected_kind)
        return validate(namespace, accession, selected_kind, source)

    monkeypatch.setattr(artifacts, "validate_file", checked)
    first = path(uri)
    second = path(uri)
    assert validations == [kind, kind]
    assert first.path == second.path
    assert first.path.read_bytes() == payload
    assert first.kind == kind
    assert first.path.is_relative_to(tmp_path)
    assert calls == [uri]
    first.path.write_bytes(b"<html>not a provider file</html>")
    with pytest.raises(ValueError):
        path(uri)


def test_failed_download_does_not_publish_cache(tmp_path, monkeypatch):
    """Publish no cache entry when a provider download fails."""
    monkeypatch.setattr(artifacts.settings, "home", tmp_path)

    def download(*args):
        args[3].write_bytes(b"partial")
        raise TimeoutError("provider stopped responding")

    monkeypatch.setattr(artifacts, "download_file", download)
    with pytest.raises(TimeoutError):
        path("pdb://1CRN")
    assert not list(tmp_path.rglob("*.cif"))


def test_short_http_body_never_becomes_a_cached_pdf(tmp_path, monkeypatch):
    """A valid magic header cannot conceal a truncated provider transfer."""
    response = Response(b"%PDF-1.7\ncut")
    response.headers["Content-Length"] = "1000"
    monkeypatch.setattr(file_sources, "urlopen", lambda *a, **kw: response)
    monkeypatch.setattr(artifacts.settings, "home", tmp_path)
    with pytest.raises(OSError, match="Incomplete HTTP response"):
        path("arxiv://2301.00001")
    assert not list(tmp_path.rglob("*.pdf"))


def test_alphafold_requires_one_model_and_validates_its_accession(
    tmp_path, monkeypatch
):
    """Accept exactly one AlphaFold model matching the requested accession."""
    payload = b"data_AF\n_ma_target_ref_db_details.db_accession P42212\n_atom_site.Cartn_x 1\n"
    responses = [
        json.dumps(
            [
                {
                    "uniprotAccession": "P42212",
                    "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P42212-F1-model_v6.cif",
                }
            ]
        ).encode(),
        payload,
    ]
    monkeypatch.setattr(
        file_sources, "urlopen", lambda *a, **kw: Response(responses.pop(0))
    )
    monkeypatch.setattr(artifacts.settings, "home", tmp_path)
    result = path("uniprot://P42212", artifact="alphafold_cif")
    assert result.path.read_bytes() == payload
    assert not responses


def _encode_metadata(content: bytes) -> dict:
    """Return the published ENCODE metadata for one file payload."""
    return {
        "accession": "ENCFF002CTW",
        "href": "/files/ENCFF002CTW/@@download/ENCFF002CTW.bed",
        "file_size": len(content),
        "md5sum": hashlib.md5(content, usedforsecurity=False).hexdigest(),
    }


def _publish_encode_file(tmp_path, monkeypatch, content: bytes = b"chr1\t0\t12\n"):
    """Publish one ENCODE file through its provider and return its cache root.

    Returns:
        Artifact cache root containing the published file and its metadata.
    """
    responses = [json.dumps(_encode_metadata(content)).encode(), content]
    monkeypatch.setattr(artifacts.settings, "home", tmp_path)
    monkeypatch.setattr(
        artifacts, "urlopen", lambda *a, **kw: Response(responses.pop(0))
    )
    published = path("encode://ENCFF002CTW")
    assert published.path.read_bytes() == content
    assert not responses
    return tmp_path / "artifacts" / "encode" / "ENCFF002CTW"


def test_encode_preserves_published_file_and_checks_checksum(tmp_path, monkeypatch):
    """Use an exact file accession, never an arbitrary experiment's first file."""
    content = b"chr1\t0\t12\n"
    responses = [json.dumps(_encode_metadata(content)).encode(), content]
    monkeypatch.setattr(artifacts.settings, "home", tmp_path)
    monkeypatch.setattr(
        artifacts, "urlopen", lambda *a, **kw: Response(responses.pop(0))
    )
    first = path("encode://ENCFF002CTW")
    assert first.path.name == "ENCFF002CTW.bed"
    assert first.path.read_bytes() == content
    assert path("encode://ENCFF002CTW").path == first.path
    with pytest.raises(FileNotFoundError, match="file accession"):
        path("encode://ENCSR163RYW")
    assert not responses


@pytest.mark.parametrize(
    "metadata",
    [
        b'{"accession":"ENCFF002CTW"}',
        b"not json",
        b'{"accession":"ENCFF002CTW","href":"/files/',
    ],
)
def test_encode_rejects_malformed_cached_metadata(metadata, tmp_path, monkeypatch):
    """Report corrupt cached metadata as one ValueError, not a lookup failure."""
    package_root = _publish_encode_file(tmp_path, monkeypatch)
    (package_root / "metadata.json").write_bytes(metadata)

    with pytest.raises(ValueError):
        path("encode://ENCFF002CTW")


@pytest.mark.parametrize(
    "published",
    [
        {"accession": "ENCFF999XXX"},
        {"href": "/files/ENCFF002CTW/@@download/../../escape.bed"},
        {"href": "https://example.org/ENCFF002CTW.bed"},
    ],
)
def test_encode_keeps_the_published_identity_check(published, tmp_path, monkeypatch):
    """Reject cached metadata that identifies another or an unsafe file."""
    package_root = _publish_encode_file(tmp_path, monkeypatch)
    metadata = {**_encode_metadata(b"chr1\t0\t12\n"), **published}
    (package_root / "metadata.json").write_text(json.dumps(metadata))

    with pytest.raises(ValueError, match="invalid file identity"):
        path("encode://ENCFF002CTW")


def test_stream_limit_applies_to_downloads_and_gzip(monkeypatch, tmp_path):
    """Stop oversized streams before publishing and retain native parser errors."""
    import gzip
    import io

    monkeypatch.setattr(file_sources.settings, "max_file_bytes", 16)
    with pytest.raises(ValueError, match="BIOV_MAX_FILE_BYTES"):
        file_sources.copy_stream(io.BytesIO(b"x" * 17), io.BytesIO())
    monkeypatch.setattr(
        file_sources, "urlopen", lambda *a, **kw: Response(gzip.compress(b"x" * 17))
    )
    with pytest.raises(ValueError, match="BIOV_MAX_FILE_BYTES"):
        file_sources.download_file(
            "emdb",
            "EMD-1",
            "density_map",
            tmp_path / "map",
            {"url": "https://example.org/map.gz"},
        )
    source = tmp_path / "record.json"
    source.write_text("{}")
    with pytest.raises(ValueError, match="Invalid study_json"):
        file_sources.validate_file(
            "clinicaltrials", "NCT00222573", "study_json", source
        )


def test_alphafold_prediction_metadata_is_bounded(monkeypatch, tmp_path):
    """Reject an oversized prediction document instead of buffering it."""
    oversized = b"[" + b" " * file_sources.MAX_PROVIDER_METADATA_BYTES
    monkeypatch.setattr(file_sources, "urlopen", lambda *a, **kw: Response(oversized))

    with pytest.raises(ValueError, match="AlphaFold prediction response is too large"):
        file_sources.download_file(
            "uniprot", "P42212", "alphafold_cif", tmp_path / "model.cif", {}
        )


@pytest.mark.parametrize(
    "bounded_document", ["PMC id conversion response", "PMC cloud metadata"]
)
def test_pmc_metadata_reads_are_bounded(bounded_document, monkeypatch, tmp_path):
    """Bound both known-small PMC metadata documents before parsing them."""
    oversized = b"{" + b" " * file_sources.MAX_PROVIDER_METADATA_BYTES
    conversion = json.dumps(
        {
            "records": [
                {
                    "pmid": 25359968,
                    "versions": [{"pmcid": "PMC4324838.1", "current": True}],
                }
            ]
        }
    ).encode()
    responses = (
        [oversized]
        if bounded_document == "PMC id conversion response"
        else [conversion, oversized]
    )
    monkeypatch.setattr(
        ncbi_files, "urlopen", lambda *a, **kw: Response(responses.pop(0))
    )

    with pytest.raises(ValueError, match=f"{bounded_document} is too large"):
        ncbi_files.download_ncbi_file(
            "pubmed", "25359968", "article_pdf", tmp_path / "article.pdf"
        )
