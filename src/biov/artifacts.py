"""Environment-local artifact paths resolved from persistent identifiers."""

import json
import os
import re
import shutil
import stat
import subprocess  # noqa: S404 - fixed argv invocation is this module's purpose
import tempfile
import zipfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import IO, Any, cast
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

from .config import settings
from .identifiers import IdentifierRef, IdentifierSyntaxError, parse_identifier

DATASETS_INCLUDE = "gff3,rna,cds,protein,genome,seq-report"
DATASETS_CLI_TIMEOUT_SECONDS = 3600.0
UNIPROT_REST_BASE_URL = "https://rest.uniprot.org"
DEFAULT_UNIPROT_TIMEOUT_SECONDS = 30.0

NCBI_DATASET_CATALOG_MEMBER = Path("ncbi_dataset/data/dataset_catalog.json")
_CATALOG_FILE_TYPE = "GENOMIC_NUCLEOTIDE_FASTA"
_MAX_CATALOG_BYTES = 8 * 1024 * 1024
_MAX_UNIPROT_FASTA_HEADER_BYTES = 64 * 1024
_MAX_UNIPROT_JSON_BYTES = 64 * 1024 * 1024
_UNIPROT_FASTA_HEADER = re.compile(rb"^>(?:sp|tr)\|([^|]+)\|")
_VERSION_SUFFIX = re.compile(r"\.[0-9]+$")
_CACHE_COMPONENT = re.compile(r"^[A-Za-z0-9._-]+$")


class ArtifactError(RuntimeError):
    """Base class for expected artifact access failures."""


class UnsupportedArtifactError(ArtifactError):
    """No provider supports the requested namespace and artifact kind."""


class ArtifactNotFoundError(ArtifactError):
    """The provider has no matching immutable biological artifact."""


class ArtifactServiceError(ArtifactError):
    """An artifact provider command or service failed."""


class ArtifactPackageError(ArtifactError):
    """A provider package violates its declared structure or integrity."""


@dataclass(frozen=True, slots=True)
class Artifact(os.PathLike[str]):
    """One original provider file exposed as a normal local path."""

    path: Path
    package_root: Path
    requested_identifier: IdentifierRef
    identifier: IdentifierRef
    kind: str
    size: int

    def __fspath__(self) -> str:
        """Return the local path for standard-library and ecosystem consumers."""
        return str(self.path)


def _run_genome_command(
    operation: str,
    accession: str,
    *arguments: str,
) -> subprocess.CompletedProcess[str]:
    """Run one official genome command without a shell.

    Returns:
        Successful completed process with captured output.

    Raises:
        ArtifactServiceError: If the CLI is missing or rejects the request.
    """
    executable = shutil.which("datasets")
    if executable is None:
        raise ArtifactServiceError(
            "install the NCBI datasets CLI; executable 'datasets' was not found on PATH"
        )
    command = [
        executable,
        operation,
        "genome",
        "accession",
        accession,
        *arguments,
    ]
    try:
        completed = subprocess.run(  # noqa: S603 - fixed argv, no shell
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=DATASETS_CLI_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as error:
        raise ArtifactServiceError(
            f"NCBI datasets CLI timed out for {accession!r}"
        ) from error
    except OSError as error:
        raise ArtifactServiceError("NCBI datasets CLI could not be executed") from error
    if completed.returncode != 0:
        diagnostic = completed.stderr.strip() or completed.stdout.strip()
        detail = f": {diagnostic}" if diagnostic else ""
        raise ArtifactServiceError(
            f"NCBI datasets CLI failed for {accession!r}{detail}"
        )
    return completed


def download_genome_package(accession: str, destination: Path) -> None:
    """Download the complete requested genome package with official CLI.

    Raises:
        ArtifactServiceError: If the CLI is missing, fails, or writes no ZIP.
    """
    _run_genome_command(
        "download",
        accession,
        "--include",
        DATASETS_INCLUDE,
        "--filename",
        str(destination),
        "--no-progressbar",
    )
    if not destination.is_file():
        raise ArtifactServiceError(
            "NCBI datasets CLI completed without writing its data package"
        )


def genome_summary(accession: str) -> str:
    """Return one validated assembled-genome JSON stdout unchanged.

    Returns:
        Native JSON text emitted by NCBI Datasets.

    Raises:
        ArtifactServiceError: If NCBI emits malformed JSON.
    """
    text = _run_genome_command("summary", accession).stdout
    try:
        summary = json.loads(text)
    except json.JSONDecodeError as error:
        raise ArtifactServiceError("NCBI genome summary is not valid JSON") from error
    if not isinstance(summary, dict):
        raise ArtifactServiceError("NCBI genome summary must be a JSON object")
    reports = summary.get("reports")
    if not isinstance(reports, list) or len(reports) != 1:
        raise ArtifactServiceError(
            "NCBI genome summary must contain exactly one report"
        )
    _matching_assembly(accession, reports)
    return text


def download_uniprot_entry(
    accession: str,
    destination: Path,
    *,
    timeout_seconds: float = DEFAULT_UNIPROT_TIMEOUT_SECONDS,
) -> None:
    """Stream an official JSON or FASTA response to its accession-named file.

    Raises:
        ValueError: If the timeout is not positive.
        ArtifactNotFoundError: If UniProt has no matching entry.
        ArtifactServiceError: If the HTTP request or local write fails.
    """
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")
    extension = destination.suffix
    media_type = {".json": "application/json", ".fasta": "text/plain"}[extension]
    request = Request(  # noqa: S310 - fixed official HTTPS origin
        f"{UNIPROT_REST_BASE_URL}/uniprotkb/{quote(accession, safe='')}{extension}",
        headers={"Accept": media_type},
    )
    try:
        with (
            urlopen(request, timeout=timeout_seconds) as response,  # noqa: S310
            destination.open("xb") as output,
        ):
            shutil.copyfileobj(response, output)
    except HTTPError as error:
        if error.code in {404, 410}:
            raise ArtifactNotFoundError(
                f"UniProt has no entry for {accession!r}"
            ) from error
        raise ArtifactServiceError(
            f"UniProt REST API failed for {accession!r}: HTTP {error.code}"
        ) from error
    except (OSError, URLError) as error:
        raise ArtifactServiceError(
            f"UniProt {extension[1:].upper()} download failed for {accession!r}"
        ) from error


def _cache_component(value: str, label: str) -> str:
    """Validate one internal cache path component.

    Returns:
        Unchanged safe component.

    Raises:
        ArtifactError: If the value is unsafe as one path component.
    """
    if _CACHE_COMPONENT.fullmatch(value) is None:
        raise ArtifactError(f"Unsafe {label} cache component {value!r}")
    return value


def _safe_package_path(value: object) -> PurePosixPath:
    """Validate a relative path from a ZIP member or dataset catalog.

    Returns:
        Safe relative package path.

    Raises:
        ArtifactPackageError: If the path is unsafe.
    """
    if not isinstance(value, str) or not value or "\\" in value:
        raise ArtifactPackageError("NCBI package contains an unsafe member path")
    parts = value.split("/")
    if (
        value.startswith("/")
        or any(part in {"", ".", ".."} for part in parts)
        or ":" in parts[0]
    ):
        raise ArtifactPackageError("NCBI package contains an unsafe member path")
    return PurePosixPath(*parts)


def _extract_complete_package(package_path: Path, destination: Path) -> None:
    """Extract every official package member without renaming or flattening.

    Raises:
        ArtifactPackageError: If the ZIP is invalid or contains unsafe paths.
    """
    try:
        with zipfile.ZipFile(package_path) as package:
            members = package.infolist()
            if not members:
                raise ArtifactPackageError("NCBI package is an empty ZIP")
            seen_paths: set[PurePosixPath] = set()
            for member in members:
                path = _safe_package_path(member.filename.removesuffix("/"))
                if stat.S_IFMT(member.external_attr >> 16) not in {
                    0,
                    stat.S_IFREG,
                    stat.S_IFDIR,
                }:
                    raise ArtifactPackageError(
                        "NCBI package contains a non-file ZIP member"
                    )
                if path in seen_paths:
                    raise ArtifactPackageError(
                        "NCBI package contains duplicate ZIP paths"
                    )
                seen_paths.add(path)
            package.extractall(destination)
    except ArtifactPackageError:
        raise
    except (OSError, RuntimeError, zipfile.BadZipFile, zipfile.LargeZipFile) as error:
        raise ArtifactPackageError(
            "NCBI genome download is not a valid ZIP package"
        ) from error


def _read_catalog(package_root: Path) -> dict[str, Any]:
    """Read the package's original dataset catalog.

    Returns:
        Parsed catalog object.

    Raises:
        ArtifactPackageError: If the catalog is absent, large, or invalid.
    """
    catalog_path = package_root / NCBI_DATASET_CATALOG_MEMBER
    try:
        if not catalog_path.is_file():
            raise ArtifactPackageError("NCBI package has no dataset catalog")
        if catalog_path.stat().st_size > _MAX_CATALOG_BYTES:
            raise ArtifactPackageError("NCBI dataset catalog is too large")
        value: Any = json.loads(catalog_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ArtifactPackageError(
            "NCBI package contains invalid dataset catalog JSON"
        ) from error
    if not isinstance(value, dict):
        raise ArtifactPackageError("NCBI dataset catalog must be an object")
    return cast("dict[str, Any]", value)


def _matching_assembly(
    requested_accession: str,
    assemblies: Any,
) -> dict[str, Any]:
    """Select the unique assembly matching a requested RefSeq accession.

    Returns:
        Matching native NCBI assembly record.

    Raises:
        ArtifactNotFoundError: If no unique requested assembly is present.
    """
    if not isinstance(assemblies, list):
        raise ArtifactNotFoundError(
            f"NCBI response does not contain exactly one assembly for "
            f"{requested_accession!r}"
        )
    accession_assemblies = [
        assembly
        for assembly in assemblies
        if isinstance(assembly, dict) and isinstance(assembly.get("accession"), str)
    ]
    if len(accession_assemblies) != 1:
        raise ArtifactNotFoundError(
            f"NCBI response does not contain exactly one assembly for "
            f"{requested_accession!r}"
        )
    assembly = accession_assemblies[0]
    accession = cast("str", assembly["accession"])
    if _VERSION_SUFFIX.search(requested_accession):
        if accession != requested_accession:
            raise ArtifactNotFoundError(
                f"NCBI did not return exact assembly version {requested_accession!r}"
            )
    elif not accession.startswith(f"{requested_accession}."):
        raise ArtifactNotFoundError(
            f"NCBI did not return a version of assembly {requested_accession!r}"
        )
    return cast("dict[str, Any]", assembly)


def _artifact_from_package(
    package_root: Path,
    requested: IdentifierRef,
    kind: str,
) -> Artifact:
    """Locate the original genomic FASTA using the official package catalog.

    Returns:
        Path-like artifact pointing at the unmodified package member.

    Raises:
        ArtifactPackageError: If package identity or content is invalid.
    """
    catalog = _read_catalog(package_root)
    assembly = _matching_assembly(requested.accession, catalog.get("assemblies"))
    accession = cast("str", assembly["accession"])
    try:
        canonical = parse_identifier(f"refseq.gcf:{accession}")
    except IdentifierSyntaxError as error:
        raise ArtifactPackageError(
            "NCBI package contains an invalid RefSeq assembly accession"
        ) from error
    files = assembly.get("files")
    if not isinstance(files, list):
        raise ArtifactPackageError("NCBI assembly catalog has no files")
    fasta_files = [
        file
        for file in files
        if isinstance(file, dict) and file.get("fileType") == _CATALOG_FILE_TYPE
    ]
    if len(fasta_files) != 1:
        raise ArtifactPackageError(
            "NCBI package must contain exactly one genomic FASTA"
        )
    fasta = fasta_files[0]
    relative_path = _safe_package_path(fasta.get("filePath"))
    if relative_path.parts[0] != canonical.accession:
        raise ArtifactPackageError(
            "NCBI catalog FASTA is outside its assembly directory"
        )
    path = package_root / "ncbi_dataset/data" / Path(*relative_path.parts)
    try:
        resolved_root = package_root.resolve(strict=True)
        resolved_path = path.resolve(strict=True)
        size = path.stat().st_size
    except OSError as error:
        raise ArtifactPackageError("NCBI catalog genomic FASTA is missing") from error
    if not resolved_path.is_relative_to(resolved_root) or not path.is_file():
        raise ArtifactPackageError("NCBI catalog genomic FASTA is unsafe")
    return Artifact(
        path=path,
        package_root=package_root,
        requested_identifier=requested,
        identifier=canonical,
        kind=kind,
        size=size,
    )


def _refseq_gcf_path(
    requested: IdentifierRef,
    kind: str,
    package_root: Path,
) -> Artifact:
    """Return the original genomic FASTA inside a complete cached package.

    Raises:
        ArtifactPackageError: If an existing cache path or package is invalid.
        ArtifactServiceError: If an atomic cache publication fails.
    """
    namespace_root = package_root.parent
    if package_root.exists():
        if not package_root.is_dir():
            raise ArtifactPackageError(
                f"artifact cache path is not a directory: {package_root}"
            )
        return _artifact_from_package(package_root, requested, kind)

    namespace_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{requested.accession}.", dir=namespace_root
    ) as staging:
        archive_path = Path(staging) / "ncbi_dataset.zip"
        extracted_root = Path(staging) / "package"
        download_genome_package(
            requested.accession,
            archive_path,
        )
        extracted_root.mkdir()
        _extract_complete_package(archive_path, extracted_root)
        _artifact_from_package(extracted_root, requested, kind)
        try:
            os.rename(extracted_root, package_root)
        except OSError as error:
            if package_root.is_dir():
                return _artifact_from_package(package_root, requested, kind)
            raise ArtifactServiceError(
                "could not publish the NCBI data package cache"
            ) from error
        return _artifact_from_package(package_root, requested, kind)


def _artifact_from_uniprot(
    path: Path,
    requested: IdentifierRef,
    kind: str,
) -> Artifact:
    """Validate one raw UniProt response and expose its original path.

    Returns:
        Path-like artifact pointing at the unchanged REST response.

    Raises:
        ArtifactNotFoundError: If the response identifies another accession.
        ArtifactPackageError: If the cached file is missing or invalid.
    """
    package_root = path.parent
    label = "entry JSON" if kind == "entry_json" else "protein FASTA"
    try:
        resolved_root = package_root.resolve(strict=True)
        resolved_path = path.resolve(strict=True)
        size = path.stat().st_size
    except OSError as error:
        raise ArtifactPackageError(f"UniProt {label} is missing") from error
    if not resolved_path.is_relative_to(resolved_root) or not path.is_file():
        raise ArtifactPackageError(f"UniProt {label} path is unsafe")
    if kind == "entry_json":
        if size > _MAX_UNIPROT_JSON_BYTES:
            raise ArtifactPackageError("UniProt entry JSON is too large")
        try:
            entry: Any = json.loads(path.read_bytes())
        except OSError as error:
            raise ArtifactPackageError("UniProt entry JSON is missing") from error
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise ArtifactPackageError("UniProt entry JSON is invalid") from error
        if not isinstance(entry, dict):
            raise ArtifactPackageError("UniProt entry JSON must be an object")
        accession = entry.get("primaryAccession")
        if accession != requested.accession:
            raise ArtifactNotFoundError(
                f"UniProt JSON did not return exact accession {requested.accession!r}"
            )
        try:
            canonical = parse_identifier(f"uniprot:{accession}")
        except IdentifierSyntaxError as error:
            raise ArtifactPackageError(
                "UniProt entry JSON contains an invalid accession"
            ) from error
    else:
        try:
            with path.open("rb") as fasta_file:
                header = fasta_file.readline(_MAX_UNIPROT_FASTA_HEADER_BYTES + 1)
        except OSError as error:
            raise ArtifactPackageError("UniProt protein FASTA is missing") from error
        if len(header) > _MAX_UNIPROT_FASTA_HEADER_BYTES:
            raise ArtifactPackageError("UniProt protein FASTA header is too large")
        match = _UNIPROT_FASTA_HEADER.match(header)
        if match is None:
            raise ArtifactPackageError("UniProt protein FASTA has an invalid header")
        try:
            accession = match.group(1).decode("ascii")
            canonical = parse_identifier(f"uniprot:{accession}")
        except (UnicodeDecodeError, IdentifierSyntaxError) as error:
            raise ArtifactPackageError(
                "UniProt protein FASTA contains an invalid accession"
            ) from error
        if accession != requested.accession:
            raise ArtifactNotFoundError(
                f"UniProt did not return exact accession {requested.accession!r}"
            )
    return Artifact(
        path=path,
        package_root=package_root,
        requested_identifier=requested,
        identifier=canonical,
        kind=kind,
        size=size,
    )


def _uniprot_path(
    requested: IdentifierRef,
    kind: str,
    package_root: Path,
) -> Artifact:
    """Return one raw JSON or FASTA response from its local cache.

    Raises:
        ArtifactPackageError: If an existing cache or response is invalid.
        ArtifactServiceError: If an atomic cache publication fails.
    """
    extension = {"entry_json": "json", "protein_fasta": "fasta"}[kind]
    filename = f"{requested.accession}.{extension}"
    try:
        package_root.mkdir(parents=True, exist_ok=True)
    except FileExistsError as error:
        raise ArtifactPackageError(
            f"artifact cache path is not a directory: {package_root}"
        ) from error
    except OSError as error:
        raise ArtifactServiceError(
            "could not create the UniProt entry cache"
        ) from error
    target = package_root / filename
    if target.exists():
        return _artifact_from_uniprot(target, requested, kind)

    with tempfile.TemporaryDirectory(
        prefix=f".{requested.accession}.{kind}.", dir=package_root.parent
    ) as staging:
        staging_path = Path(staging) / filename
        download_uniprot_entry(requested.accession, staging_path)
        _artifact_from_uniprot(staging_path, requested, kind)
        try:
            os.replace(staging_path, target)
        except OSError as error:
            if target.is_file():
                return _artifact_from_uniprot(target, requested, kind)
            raise ArtifactServiceError(
                f"could not publish the UniProt {kind} cache"
            ) from error
        return _artifact_from_uniprot(target, requested, kind)


ARTIFACT_PROVIDERS: dict[
    tuple[str, str], Callable[[IdentifierRef, str, Path], Artifact]
] = {
    ("refseq.gcf", "genome_fasta"): _refseq_gcf_path,
    ("uniprot", "protein_fasta"): _uniprot_path,
    ("uniprot", "entry_json"): _uniprot_path,
}

DEFAULT_ARTIFACT_KINDS = {
    "refseq.gcf": "genome_fasta",
    "uniprot": "protein_fasta",
}


def path(
    identifier: str | IdentifierRef,
    *,
    artifact: str | None = None,
) -> Artifact:
    """Return an identifier-backed artifact path in the current environment.

    The returned object is a normal os.PathLike and can be passed directly
    to Biopython and other libraries. Concrete paths intentionally remain local
    to the process and executor running this function.

    Args:
        identifier: Exact identifier reference or an already parsed reference.
        artifact: Provider-independent artifact kind. If omitted, use the
            namespace's default kind.

    Returns:
        Original file inside the complete cached provider package.

    Raises:
        UnsupportedArtifactError: If no provider supports the pair.
    """
    requested = (
        identifier
        if isinstance(identifier, IdentifierRef)
        else parse_identifier(identifier)
    )
    kind = (
        artifact
        if artifact is not None
        else DEFAULT_ARTIFACT_KINDS.get(requested.namespace, "<default>")
    )
    provider = ARTIFACT_PROVIDERS.get((requested.namespace, kind))
    if provider is None:
        raise UnsupportedArtifactError(
            f"No artifact provider for namespace {requested.namespace!r} "
            f"and artifact {kind!r}"
        )
    package_root = (
        settings.home
        / "artifacts"
        / _cache_component(requested.namespace, "namespace")
        / _cache_component(requested.accession, "accession")
    )
    return provider(requested, kind, package_root)


def open(
    identifier: str | IdentifierRef,
    *,
    artifact: str | None = None,
    mode: str = "rb",
    encoding: str | None = None,
) -> IO[Any]:
    """Open an identifier-backed artifact through the read-only file API.

    Args:
        identifier: Exact identifier reference or parsed reference.
        artifact: Provider-independent artifact kind. If omitted, use the
            namespace's default kind.
        mode: Read-only text or binary file mode.
        encoding: Text encoding; defaults to UTF-8 for text mode.

    Returns:
        A standard file object for the environment-local artifact.

    Raises:
        ValueError: If a mutating or unsupported mode is requested.
    """
    if not mode.startswith("r") or any(flag in mode for flag in "wax+"):
        raise ValueError("open supports read-only modes")
    if encoding is None and "b" not in mode:
        encoding = "utf-8"
    return path(identifier, artifact=artifact).path.open(mode, encoding=encoding)


__all__ = [
    "DATASETS_CLI_TIMEOUT_SECONDS",
    "DATASETS_INCLUDE",
    "DEFAULT_UNIPROT_TIMEOUT_SECONDS",
    "NCBI_DATASET_CATALOG_MEMBER",
    "UNIPROT_REST_BASE_URL",
    "Artifact",
    "ArtifactError",
    "ArtifactNotFoundError",
    "ArtifactPackageError",
    "ArtifactServiceError",
    "UnsupportedArtifactError",
    "open",
    "path",
]
