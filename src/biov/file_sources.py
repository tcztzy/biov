"""Download documented sequence, structure, article, and study files."""

import gzip
import json
import re
from http.client import HTTPResponse
from importlib.metadata import version
from pathlib import Path
from typing import Any, BinaryIO
from urllib.parse import quote, urlsplit
from urllib.request import Request, urlopen
from xml.etree import ElementTree  # noqa: S405 - no external entity resolution

from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from .config import settings

# Provider metadata documents are known to be small; artifact bodies stay
# bounded only by the operator's optional per-file limit.
MAX_PROVIDER_METADATA_BYTES = 1024 * 1024


def copy_stream(source: BinaryIO | gzip.GzipFile, output: BinaryIO) -> None:
    """Copy bytes with the operator's optional per-file size limit.

    Raises:
        ValueError: If the stream exceeds the configured size limit.
    """
    total = 0
    while chunk := source.read(1024 * 1024):
        total += len(chunk)
        if settings.max_file_bytes is not None and total > settings.max_file_bytes:
            raise ValueError("File exceeds BIOV_MAX_FILE_BYTES")
        output.write(chunk)


def copy_response(response: HTTPResponse, output: BinaryIO) -> None:
    """Copy a provider response and reject an incomplete declared HTTP body.

    Args:
        response: Open provider response with its HTTP headers.
        output: Seekable staging file receiving the unchanged bytes.

    Raises:
        OSError: If the received length differs from Content-Length.
    """
    expected = response.headers.get("Content-Length")
    start = output.tell()
    copy_stream(response, output)
    received = output.tell() - start
    if expected is not None and received != int(expected):
        raise OSError(
            f"Incomplete HTTP response: received {received} of {expected} bytes"
        )


def load_metadata_json(response: HTTPResponse, label: str) -> Any:
    """Read one known-small provider metadata document as JSON.

    Args:
        response: Open provider response carrying the document.
        label: Document name used in the size error.

    Returns:
        Parsed provider JSON document.

    Raises:
        ValueError: If the document exceeds the fixed metadata bound.
    """
    raw = response.read(MAX_PROVIDER_METADATA_BYTES + 1)
    if len(raw) > MAX_PROVIDER_METADATA_BYTES:
        raise ValueError(f"{label} is too large")
    return json.loads(raw)


def download_file(
    namespace: str, accession: str, kind: str, destination: Path, properties: dict
) -> None:
    """Download one declared file representation, decompressing transport gzip.

    Args:
        namespace: Validated identifiers.org namespace.
        accession: Validated namespace-local accession.
        kind: File kind from the packaged capability manifest.
        destination: New staging file owned by the artifact cache.
        properties: Packaged file URL, suffix, and media type.

    Raises:
        FileNotFoundError: If AlphaFold has no unique matching model.
        ValueError: If the model URL is outside the official file host, or the
            prediction metadata exceeds the fixed metadata bound.
        OSError: If the response body is truncated or exceeds the size limit.
    """
    if kind == "alphafold_cif":
        with urlopen(
            f"https://alphafold.ebi.ac.uk/api/prediction/{quote(accession, safe='')}",
            timeout=60,
        ) as response:
            models = load_metadata_json(response, "AlphaFold prediction response")
        if len(models) != 1 or models[0]["uniprotAccession"] != accession:
            raise FileNotFoundError(
                f"AlphaFold has no unique model for {accession!r}; select a model explicitly"
            )
        url = models[0]["cifUrl"]
        parsed = urlsplit(url)
        if (
            parsed.scheme != "https"
            or parsed.netloc != "alphafold.ebi.ac.uk"
            or not parsed.path.startswith("/files/")
            or not parsed.path.endswith(".cif")
        ):
            raise ValueError("Unexpected AlphaFold coordinate URL")
    else:
        url = properties["url"].format(
            accession=quote(accession, safe="/" if namespace == "arxiv" else ""),
            number=accession.removeprefix("EMD-"),
        )
    request = Request(url, headers={"User-Agent": f"BioV/{version('biov')}"})  # noqa: S310
    with urlopen(request, timeout=60) as response, destination.open("xb") as output:  # noqa: S310
        if url.endswith(".gz"):
            with gzip.GzipFile(fileobj=response) as content:
                try:
                    copy_stream(content, output)
                except EOFError as error:
                    raise OSError("Truncated gzip response body") from error
        else:
            copy_response(response, output)


def validate_file(namespace: str, accession: str, kind: str, source: Path) -> None:
    """Check file format and provider identity where included in the format.

    Args:
        namespace: Validated identifiers.org namespace.
        accession: Requested accession.
        kind: Declared file representation.
        source: Downloaded or cached uncompressed file.

    Raises:
        ValueError: If the parsed file fails format or identity checks.

    Native parser and I/O errors propagate unchanged.
    """
    valid = False
    if kind in {"structure_cif", "alphafold_cif"}:
        record = MMCIF2Dict(str(source))
        if kind == "structure_cif":
            valid = record.get("_entry.id", [""])[0].upper() == accession.upper()
        else:
            valid = accession in record.get(
                "_ma_target_ref_db_details.db_accession", []
            )
        valid = valid and bool(record.get("_atom_site.Cartn_x"))
    elif kind == "sequence_fasta":
        record = SeqIO.read(source, "fasta")
        valid = bool(record.seq) and (
            record.id == accession
            or ("." not in accession and record.id.rsplit(".", 1)[0] == accession)
        )
    elif kind == "structure_sdf":
        text = source.read_text()
        marker = (
            rf">\s*<chembl_id>\s*\n{re.escape(accession)}\s*\n"
            if namespace == "chembl.compound"
            else rf"\A{re.escape(accession)}\s*\n"
        )
        valid = (
            re.search(marker, text) is not None
            and "M  END" in text
            # ChEMBL's single-record response ends after its SDF properties.
            and (namespace == "chembl.compound" or text.rstrip().endswith("$$$$"))
        )
    elif kind == "study_json":
        with source.open() as content:
            record = json.load(content)
        section = record.get("protocolSection") if isinstance(record, dict) else None
        module = (
            section.get("identificationModule") if isinstance(section, dict) else None
        )
        valid = isinstance(module, dict) and module.get("nctId") == accession
    elif kind in {"label_xml", "pathway_sbml"}:
        root = ElementTree.parse(source).getroot()  # noqa: S314
        if kind == "label_xml":
            identity = root.find("{urn:hl7-org:v3}setId")
            valid = identity is not None and identity.get("root") == accession
        else:
            valid = root.tag.endswith("}sbml") and any(
                accession in value
                for node in root.iter()
                for value in node.attrib.values()
            )
    else:
        with source.open("rb") as content:
            header = content.read(1024)
        if kind == "article_pdf":
            valid = header.startswith(b"%PDF-")
        elif kind == "density_map":
            valid = len(header) == 1024 and header[208:212] == b"MAP "
    if not valid:
        raise ValueError(f"Invalid {kind} file for {namespace}:{accession}")
