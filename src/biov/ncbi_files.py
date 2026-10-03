"""Download native NCBI records, GEO tables, and available PMC article PDFs."""

import gzip
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import parse_qs, urlencode, urlsplit
from urllib.request import urlopen
from xml.etree import ElementTree  # noqa: S405 - no external entity resolution

import fsspec

from .file_sources import copy_response, copy_stream, load_metadata_json


def _pubmed_pdf(accession: str, destination: Path) -> None:
    """Resolve the current PMC article version and download its licensed PDF.

    Raises:
        FileNotFoundError: If no current PMC version or PDF exists.
        ValueError: If the PDF identity or checksum is invalid, or a metadata
            document exceeds the fixed metadata bound.
    """
    query = urlencode(
        {
            "ids": accession,
            "idtype": "pmid",
            "format": "json",
            "versions": "yes",
            "tool": "biov",
        }
    )
    with urlopen(
        f"https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/?{query}",
        timeout=60,
    ) as response:
        records = load_metadata_json(response, "PMC id conversion response")["records"]
    if len(records) != 1 or str(records[0].get("pmid")) != accession:
        raise FileNotFoundError(f"PubMed {accession} has no matching PMC article")
    versions = [
        version["pmcid"]
        for version in records[0].get("versions", [])
        if version.get("current") is True
    ]
    if len(versions) != 1 or not re.fullmatch(r"PMC\d+\.\d+", versions[0]):
        raise FileNotFoundError(f"PubMed {accession} has no unique current PMC version")
    version = versions[0]
    with urlopen(
        f"https://pmc-oa-opendata.s3.amazonaws.com/metadata/{version}.json", timeout=60
    ) as response:
        metadata = load_metadata_json(response, "PMC cloud metadata")
    if (
        str(metadata["pmid"]) != accession
        or f"{metadata['pmcid']}.{metadata['version']}" != version
    ):
        raise ValueError("PMC metadata does not match the requested article")
    if not metadata.get("pdf_url"):
        raise FileNotFoundError(f"PubMed {accession} has no downloadable PMC PDF")
    pdf = urlsplit(metadata["pdf_url"])
    if (
        pdf.scheme != "s3"
        or pdf.netloc != "pmc-oa-opendata"
        or pdf.path != f"/{version}/{version}.pdf"
    ):
        raise ValueError("PMC metadata contains an unexpected PDF URL")
    checksum = parse_qs(pdf.query)["md5"]
    with (
        urlopen(
            f"https://pmc-oa-opendata.s3.amazonaws.com{pdf.path}", timeout=60
        ) as response,
        destination.open("xb") as output,
    ):
        copy_response(response, output)
    with destination.open("rb") as content:
        actual = hashlib.file_digest(
            content, lambda: hashlib.md5(usedforsecurity=False)
        ).hexdigest()
    if checksum != [actual]:
        raise ValueError("PMC PDF does not match its published MD5 checksum")


def _geo_file(accession: str, kind: str, destination: Path) -> None:
    """Download one GEO Series matrix or an accession's complete SOFT record.

    Raises:
        FileNotFoundError: If the accession has no unique Series matrix.
    """
    prefix = accession[:3]
    group = f"{accession[:-3] if len(accession) > 6 else prefix}nnn"
    if kind == "expression_matrix":
        if prefix != "GSE":
            raise FileNotFoundError(
                f"{accession} is not a GEO Series; request artifact='soft' for its data"
            )
        directory = (
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/{group}/{accession}/matrix/"
        )
        candidates = [
            url
            for url in fsspec.filesystem("https").ls(directory, detail=False)
            if re.fullmatch(
                re.escape(directory + accession)
                + r"(?:-GPL\d+)?_series_matrix\.txt\.gz",
                url,
            )
        ]
        if len(candidates) != 1:
            raise FileNotFoundError(
                f"{accession} has {len(candidates)} Series matrices; "
                f"choose a file at {directory}"
            )
        url = candidates[0]
    elif prefix in {"GSE", "GDS"}:
        directory = "series" if prefix == "GSE" else "datasets"
        suffix = "_family.soft.gz" if prefix == "GSE" else "_full.soft.gz"
        url = (
            f"https://ftp.ncbi.nlm.nih.gov/geo/{directory}/{group}/"
            f"{accession}/soft/{accession}{suffix}"
        )
    else:
        query = urlencode(
            {"acc": accession, "targ": "self", "view": "full", "form": "text"}
        )
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?{query}"
    with (
        urlopen(url, timeout=60) as response,  # noqa: S310 - fixed NCBI URLs
        destination.open("xb") as output,
    ):
        if url.endswith(".gz"):
            with gzip.GzipFile(fileobj=response) as content:
                copy_stream(content, output)
        else:
            copy_response(response, output)


def validate_ncbi_file(namespace: str, accession: str, kind: str, source: Path) -> None:
    """Check the downloaded representation and its accession before cache reuse.

    Args:
        namespace: Validated identifier namespace.
        accession: Validated accession in that namespace.
        kind: Artifact kind selected from the capability manifest.
        source: Local uncompressed data file.

    Raises:
        ValueError: If parsed data is mismatched or lacks its table.

    Native parser and I/O errors propagate unchanged.
    """
    valid = False
    if namespace == "pubmed":
        with source.open("rb") as content:
            valid = content.read(5) == b"%PDF-"
    elif namespace == "dbsnp":
        with source.open("rb") as content:
            record = json.load(content)
        valid = isinstance(record, dict) and (
            str(record.get("refsnp_id")) == accession.removeprefix("rs")
        )
    elif namespace == "clinvar":
        # ElementTree does not resolve external entities; malformed XML fails.
        records = ElementTree.parse(source).getroot().findall("VariationArchive")  # noqa: S314
        valid = len(records) == 1 and records[0].get("VariationID") == accession
    else:
        marker = (
            "!series_matrix_table_begin"
            if kind == "expression_matrix"
            else {
                "GSE": "^SERIES = ",
                "GDS": "^DATASET = ",
                "GPL": "^PLATFORM = ",
                "GSM": "^SAMPLE = ",
            }[accession[:3]]
            + accession
        )
        with source.open(encoding="utf-8") as content:
            matching_accession = kind != "expression_matrix"
            for line in content:
                if line.rstrip("\r\n") == f'!Series_geo_accession\t"{accession}"':
                    matching_accession = True
                if line.rstrip("\r\n") == marker:
                    valid = matching_accession
                    break
            if valid and kind == "expression_matrix":
                # Empty RNA-seq matrices only carry metadata, not expression data.
                header = next(content, "")
                row = next(content, "")
                valid = (
                    header.startswith('"ID_REF"\t')
                    and bool(row.strip())
                    and not row.startswith("!")
                )
    if not valid:
        raise ValueError(
            f"{namespace}:{accession} did not return a matching {kind} file"
        )


def download_ncbi_file(
    namespace: str, accession: str, kind: str, destination: Path
) -> None:
    """Download and validate one NCBI data file without replacing another format.

    Args:
        namespace: Namespace already validated by the identifier parser.
        accession: Accession already validated with the identifiers.org registry.
        kind: Supported artifact kind selected by the capability manifest.
        destination: Nonexistent path in the caller's atomic staging directory.

    """
    if namespace == "pubmed":
        _pubmed_pdf(accession, destination)
    elif namespace == "geo":
        _geo_file(accession, kind, destination)
    else:
        if namespace == "clinvar":
            query = urlencode(
                {
                    "db": "clinvar",
                    "id": accession,
                    "rettype": "vcv",
                    "is_variationid": "true",
                    "retmode": "xml",
                }
            )
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}"
        else:
            url = (
                "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/"
                + accession.removeprefix("rs")
            )
        with (
            urlopen(url, timeout=60) as response,  # noqa: S310 - fixed official URLs
            destination.open("xb") as output,
        ):
            copy_response(response, output)
