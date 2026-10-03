"""Environment-local artifact paths resolved from persistent identifiers."""

import hashlib
import json
import os
import re
import shutil
import stat
import subprocess  # noqa: S404 - fixed argv invocation is this module's purpose
import tempfile
import zipfile
from collections.abc import Callable
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import IO, Any, cast
from urllib.parse import quote
from urllib.request import Request, urlopen

from pydantic import BaseModel

from .capabilities import load_capabilities as _artifact_capabilities
from .config import settings
from .file_sources import copy_response, download_file, validate_file
from .identifiers import IdentifierRef, IdentifierSyntaxError, parse_identifier
from .ncbi_files import download_ncbi_file, validate_ncbi_file

DATASETS_INCLUDE = "gff3,rna,cds,protein,genome,seq-report"
DATASETS_CLI_TIMEOUT_SECONDS = 3600.0
UNIPROT_REST_BASE_URL = "https://rest.uniprot.org"
DEFAULT_UNIPROT_TIMEOUT_SECONDS = 30.0

NCBI_DATASET_CATALOG_MEMBER = Path("ncbi_dataset/data/dataset_catalog.json")
_MAX_CATALOG_BYTES = 8 * 1024 * 1024
_MAX_UNIPROT_FASTA_HEADER_BYTES = 64 * 1024
_MAX_UNIPROT_JSON_BYTES = 64 * 1024 * 1024
_UNIPROT_FASTA_HEADER = re.compile(rb"^>(?:sp|tr)\|([^|]+)\|")
_VERSION_SUFFIX = re.compile(r"\.[0-9]+$")


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
        RuntimeError: If the CLI is missing or rejects the request.
    """
    executable = shutil.which("datasets")
    if executable is None:
        raise RuntimeError(
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
    completed = subprocess.run(  # noqa: S603 - fixed argv, no shell
        command,
        check=False,
        capture_output=True,
        text=True,
        timeout=DATASETS_CLI_TIMEOUT_SECONDS,
    )
    if completed.returncode != 0:
        diagnostic = completed.stderr.strip() or completed.stdout.strip()
        detail = f": {diagnostic}" if diagnostic else ""
        raise RuntimeError(f"NCBI datasets CLI failed for {accession!r}{detail}")
    return completed


def download_genome_package(accession: str, destination: Path) -> None:
    """Download the complete requested genome package with official CLI.

    Raises:
        RuntimeError: If the CLI is missing, fails, or writes no ZIP.
        ValueError: If the package exceeds ``BIOV_MAX_FILE_BYTES``.
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
        raise RuntimeError(
            "NCBI datasets CLI completed without writing its data package"
        )
    limit = settings.max_file_bytes
    if limit is not None and destination.stat().st_size > limit:
        raise ValueError("NCBI Datasets package exceeds BIOV_MAX_FILE_BYTES")


def genome_summary(accession: str) -> str:
    """Return one validated assembled-genome JSON stdout unchanged.

    Returns:
        Native JSON text emitted by NCBI Datasets.

    Raises:
        ValueError: If NCBI emits malformed JSON or a non-object summary.
    """
    text = _run_genome_command("summary", accession).stdout
    summary = json.loads(text)
    if not isinstance(summary, dict):
        raise ValueError(  # noqa: TRY004 - untrusted provider JSON, not a caller type error.
            "NCBI genome summary must be a JSON object"
        )
    reports = summary.get("reports")
    if not isinstance(reports, list) or len(reports) != 1:
        raise ValueError("NCBI genome summary must contain exactly one report")
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
    """
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")
    extension = destination.suffix
    media_type = {".json": "application/json", ".fasta": "text/plain"}[extension]
    request = Request(  # noqa: S310 - fixed official HTTPS origin
        f"{UNIPROT_REST_BASE_URL}/uniprotkb/{quote(accession, safe='')}{extension}",
        headers={"Accept": media_type},
    )
    with (
        urlopen(request, timeout=timeout_seconds) as response,  # noqa: S310
        destination.open("xb") as output,
    ):
        copy_response(response, output)


def _cache_component(value: str, label: str) -> str:
    """Validate one internal cache path component.

    Returns:
        Percent-encoded safe component.

    Raises:
        ValueError: If the value is unsafe as one path component.
    """
    if not value or value in {".", ".."}:
        raise ValueError(f"Unsafe {label} cache component {value!r}")
    return quote(value, safe="")


def _safe_package_path(value: object) -> PurePosixPath:
    """Validate a relative path from a ZIP member or dataset catalog.

    Returns:
        Safe relative package path.

    Raises:
        ValueError: If the path is unsafe.
    """
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError("NCBI package contains an unsafe member path")
    parts = value.split("/")
    if (
        value.startswith("/")
        or any(part in {"", ".", ".."} for part in parts)
        or ":" in parts[0]
    ):
        raise ValueError("NCBI package contains an unsafe member path")
    return PurePosixPath(*parts)


def _extract_complete_package(package_path: Path, destination: Path) -> None:
    """Extract every official package member without renaming or flattening.

    Each member is checked against the operator's per-file ceiling before the
    archive is unpacked, so an oversized decompressed file is rejected instead
    of published.

    Raises:
        ValueError: If the ZIP contains unsafe or duplicate member paths, or a
            member exceeds ``BIOV_MAX_FILE_BYTES``.
    """
    limit = settings.max_file_bytes
    with zipfile.ZipFile(package_path) as package:
        members = package.infolist()
        if not members:
            raise ValueError("NCBI package is an empty ZIP")
        seen_paths: set[PurePosixPath] = set()
        for member in members:
            path = _safe_package_path(member.filename.removesuffix("/"))
            if stat.S_IFMT(member.external_attr >> 16) not in {
                0,
                stat.S_IFREG,
                stat.S_IFDIR,
            }:
                raise ValueError("NCBI package contains a non-file ZIP member")
            if path in seen_paths:
                raise ValueError("NCBI package contains duplicate ZIP paths")
            if limit is not None and member.file_size > limit:
                raise ValueError(
                    f"NCBI package member {path} exceeds BIOV_MAX_FILE_BYTES"
                )
            seen_paths.add(path)
        package.extractall(destination)


def _read_catalog(package_root: Path) -> dict[str, Any]:
    """Read the package's original dataset catalog.

    Returns:
        Parsed catalog object.

    Raises:
        ValueError: If the catalog is absent, large, invalid, or not an object.
    """
    catalog_path = package_root / NCBI_DATASET_CATALOG_MEMBER
    if not catalog_path.is_file():
        raise ValueError("NCBI package has no dataset catalog")
    if catalog_path.stat().st_size > _MAX_CATALOG_BYTES:
        raise ValueError("NCBI dataset catalog is too large")
    value: Any = json.loads(catalog_path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(  # noqa: TRY004 - untrusted provider JSON, not a caller type error.
            "NCBI dataset catalog must be an object"
        )
    return cast("dict[str, Any]", value)


def _matching_assembly(
    requested_accession: str,
    assemblies: Any,
) -> dict[str, Any]:
    """Select the unique assembly matching a requested RefSeq accession.

    Returns:
        Matching native NCBI assembly record.

    Raises:
        FileNotFoundError: If no unique requested assembly is present.
    """
    if not isinstance(assemblies, list):
        raise FileNotFoundError(
            f"NCBI response does not contain exactly one assembly for "
            f"{requested_accession!r}"
        )
    accession_assemblies = [
        assembly
        for assembly in assemblies
        if isinstance(assembly, dict) and isinstance(assembly.get("accession"), str)
    ]
    if len(accession_assemblies) != 1:
        raise FileNotFoundError(
            f"NCBI response does not contain exactly one assembly for "
            f"{requested_accession!r}"
        )
    assembly = accession_assemblies[0]
    accession = cast("str", assembly["accession"])
    if _VERSION_SUFFIX.search(requested_accession):
        if accession != requested_accession:
            raise FileNotFoundError(
                f"NCBI did not return exact assembly version {requested_accession!r}"
            )
    elif not accession.startswith(f"{requested_accession}."):
        raise FileNotFoundError(
            f"NCBI did not return a version of assembly {requested_accession!r}"
        )
    return cast("dict[str, Any]", assembly)


@dataclass(frozen=True, slots=True)
class _CatalogMember:
    """One catalog-selected member of a complete NCBI data package."""

    relative_path: PurePosixPath
    identifier: IdentifierRef
    label: str


def _resolve_package_member(
    package_root: Path,
    requested: IdentifierRef,
    kind: str,
) -> _CatalogMember:
    """Read one package catalog and select the unique requested member.

    Returns:
        Catalog-resolved member path and canonical package identity.

    Raises:
        ValueError: If the catalog lacks the requested unique member.
    """
    catalog = _read_catalog(package_root)
    assembly = _matching_assembly(requested.accession, catalog.get("assemblies"))
    accession = cast("str", assembly["accession"])
    try:
        canonical = parse_identifier(f"refseq.gcf:{accession}")
    except IdentifierSyntaxError as error:
        raise ValueError(
            "NCBI package contains an invalid RefSeq assembly accession"
        ) from error
    capability = _artifact_capabilities()["namespaces"]["refseq.gcf"]["kinds"][kind]
    file_type, label = capability["fileType"], capability["label"]
    files = assembly.get("files")
    if not isinstance(files, list):
        raise ValueError(  # noqa: TRY004 - untrusted provider JSON, not a caller type error.
            "NCBI assembly catalog has no files"
        )
    matching_files = [
        file
        for file in files
        if isinstance(file, dict) and file.get("fileType") == file_type
    ]
    if len(matching_files) != 1:
        raise ValueError(f"NCBI package must contain exactly one {label}")
    relative_path = _safe_package_path(matching_files[0].get("filePath"))
    if relative_path.parts[0] != canonical.accession:
        raise ValueError(f"NCBI catalog {label} is outside its assembly directory")
    return _CatalogMember(relative_path, canonical, label)


def _artifact_from_root(
    package_root: Path,
    requested: IdentifierRef,
    kind: str,
    member: _CatalogMember,
) -> Artifact:
    """Validate one resolved catalog member under a package root.

    Returns:
        Path-like artifact pointing at the unmodified package member.

    Raises:
        ValueError: If the member is missing or escapes its package root.
    """
    label = member.label
    path = package_root / "ncbi_dataset/data" / Path(*member.relative_path.parts)
    try:
        resolved_root = package_root.resolve(strict=True)
        resolved_path = path.resolve(strict=True)
        size = path.stat().st_size
    except OSError as error:
        raise ValueError(f"NCBI catalog {label} is missing") from error
    if not resolved_path.is_relative_to(resolved_root) or not path.is_file():
        raise ValueError(f"NCBI catalog {label} is unsafe")
    return Artifact(
        path=path,
        package_root=package_root,
        requested_identifier=requested,
        identifier=member.identifier,
        kind=kind,
        size=size,
    )


def _refseq_gcf_path(
    requested: IdentifierRef,
    kind: str,
    package_root: Path,
) -> Artifact:
    """Return the original requested member inside a complete cached package.

    The catalog is read once per request: the resolved member is validated
    against the staging root and then against the published root.

    Raises:
        ValueError: If an existing cache path or package is invalid.
        OSError: If an atomic cache publication fails.
    """
    namespace_root = package_root.parent
    if package_root.exists():
        if not package_root.is_dir():
            raise ValueError(f"artifact cache path is not a directory: {package_root}")
        member = _resolve_package_member(package_root, requested, kind)
        return _artifact_from_root(package_root, requested, kind, member)

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
        member = _resolve_package_member(extracted_root, requested, kind)
        # The staged member must exist before the package is published.
        _artifact_from_root(extracted_root, requested, kind, member)
        try:
            os.rename(extracted_root, package_root)
        except OSError:
            if package_root.is_dir():
                # Another request published its own package first.
                published = _resolve_package_member(package_root, requested, kind)
                return _artifact_from_root(package_root, requested, kind, published)
            raise
        return _artifact_from_root(package_root, requested, kind, member)


def _artifact_from_uniprot(
    path: Path,
    requested: IdentifierRef,
    kind: str,
) -> Artifact:
    """Validate one raw UniProt response and expose its original path.

    Returns:
        Path-like artifact pointing at the unchanged REST response.

    Raises:
        FileNotFoundError: If the file is missing or identifies another accession.
        ValueError: If the cached content or path is invalid.
    """
    package_root = path.parent
    label = "entry JSON" if kind == "entry_json" else "protein FASTA"
    resolved_root = package_root.resolve(strict=True)
    resolved_path = path.resolve(strict=True)
    size = path.stat().st_size
    if not resolved_path.is_relative_to(resolved_root) or not path.is_file():
        raise ValueError(f"UniProt {label} path is unsafe")
    if kind == "entry_json":
        if size > _MAX_UNIPROT_JSON_BYTES:
            raise ValueError("UniProt entry JSON is too large")
        entry: Any = json.loads(path.read_bytes())
        if not isinstance(entry, dict):
            raise ValueError("UniProt entry JSON must be an object")
        accession = entry.get("primaryAccession")
        if accession != requested.accession:
            raise FileNotFoundError(
                f"UniProt JSON did not return exact accession {requested.accession!r}"
            )
        canonical = parse_identifier(f"uniprot:{accession}")
    else:
        with path.open("rb") as fasta_file:
            header = fasta_file.readline(_MAX_UNIPROT_FASTA_HEADER_BYTES + 1)
        if len(header) > _MAX_UNIPROT_FASTA_HEADER_BYTES:
            raise ValueError("UniProt protein FASTA header is too large")
        match = _UNIPROT_FASTA_HEADER.match(header)
        if match is None:
            raise ValueError("UniProt protein FASTA has an invalid header")
        accession = match.group(1).decode("ascii")
        canonical = parse_identifier(f"uniprot:{accession}")
        if accession != requested.accession:
            raise FileNotFoundError(
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
        OSError: If an atomic cache publication fails.
    """
    if kind == "alphafold_cif":
        return _file_path(requested, kind, package_root)
    extension = {"entry_json": "json", "protein_fasta": "fasta"}[kind]
    filename = f"{requested.accession}.{extension}"
    package_root.mkdir(parents=True, exist_ok=True)
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
        except OSError:
            if target.is_file():
                return _artifact_from_uniprot(target, requested, kind)
            raise
        return _artifact_from_uniprot(target, requested, kind)


def _file_path(requested: IdentifierRef, kind: str, package_root: Path) -> Artifact:
    """Download, validate, and atomically cache one declared analysis file.

    Returns:
        Path-like artifact for the validated file in the current environment.

    Raises:
        ValueError: If the cache contains an unsafe or empty file.
    """
    record = _artifact_capabilities()["namespaces"][requested.namespace]
    properties = record["kinds"][kind]
    target = (
        package_root
        / f"{_cache_component(requested.accession, 'accession')}{properties['suffix']}"
    )
    validate = (
        validate_ncbi_file if record["provider"] == "ncbi_file" else validate_file
    )
    package_root.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        with tempfile.TemporaryDirectory(
            prefix=".download-", dir=package_root
        ) as staging:
            staged = Path(staging) / target.name
            if record["provider"] == "ncbi_file":
                download_ncbi_file(
                    requested.namespace, requested.accession, kind, staged
                )
            else:
                download_file(
                    requested.namespace, requested.accession, kind, staged, properties
                )
            validate(requested.namespace, requested.accession, kind, staged)
            os.replace(staged, target)
    else:
        validate(requested.namespace, requested.accession, kind, target)
    info = target.lstat()
    if not stat.S_ISREG(info.st_mode) or info.st_size == 0:
        raise ValueError(f"Invalid cached file: {target}")
    return Artifact(target, package_root, requested, requested, kind, info.st_size)


class _EncodeFileMetadata(BaseModel):
    """Published ENCODE file fields that identify and verify one download."""

    accession: str
    href: str
    file_size: int
    md5sum: str


def _encode_href(metadata: _EncodeFileMetadata, accession: str) -> str:
    """Validate the provider identity before constructing paths or URLs.

    Returns:
        Validated relative download URL.

    Raises:
        ValueError: If the URL does not identify the requested file.
    """
    prefix = f"/files/{accession}/@@download/{accession}."
    href = metadata.href
    if (
        metadata.accession != accession
        or re.fullmatch(re.escape(prefix) + r"[A-Za-z0-9.]+", href) is None
    ):
        raise ValueError("ENCODE metadata has an invalid file identity")
    return href


def _encode_artifact(package_root: Path, requested: IdentifierRef) -> Artifact:
    """Check one ENCODE file against its original download metadata.

    Returns:
        Original file with the extension provided by ENCODE.

    Raises:
        ValueError: If metadata identifies an unsafe or incomplete file.
    """
    metadata = _EncodeFileMetadata.model_validate_json(
        (package_root / "metadata.json").read_bytes()
    )
    href = _encode_href(metadata, requested.accession)
    target = package_root / PurePosixPath(href).name
    if (
        target.is_symlink()
        or not target.is_file()
        or target.stat().st_size != metadata.file_size
    ):
        raise ValueError("ENCODE file does not match its published size")
    return Artifact(
        target, package_root, requested, requested, "data_file", metadata.file_size
    )


def _encode_path(requested: IdentifierRef, kind: str, package_root: Path) -> Artifact:
    """Cache a specific ENCODE file and verify its published MD5 on download.

    Returns:
        Original data file, preserving its name, format, and compression.

    Raises:
        FileNotFoundError: If an experiment or other non-file ID is supplied.
        ValueError: If metadata is oversized or the file fails its checksum.
        OSError: If cache publication fails.
    """
    if not requested.accession.startswith("ENCFF"):
        raise FileNotFoundError(
            "Choose an ENCODE file accession (ENCFF); an experiment has multiple files"
        )
    if package_root.exists():
        return _encode_artifact(package_root, requested)
    package_root.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".encode-", dir=package_root.parent
    ) as staging:
        staged = Path(staging) / "package"
        staged.mkdir()
        with urlopen(
            f"https://www.encodeproject.org/files/{requested.accession}/?format=json",
            timeout=60,
        ) as response:
            raw = response.read(_MAX_CATALOG_BYTES + 1)
        if len(raw) > _MAX_CATALOG_BYTES:
            raise ValueError("ENCODE metadata is too large")
        metadata = _EncodeFileMetadata.model_validate_json(raw)
        href = _encode_href(metadata, requested.accession)
        (staged / "metadata.json").write_bytes(raw)
        target = staged / PurePosixPath(href).name
        with (
            urlopen("https://www.encodeproject.org" + href, timeout=60) as response,
            target.open("xb") as output,
        ):
            copy_response(response, output)
        with target.open("rb") as content:
            digest = hashlib.file_digest(
                content, lambda: hashlib.md5(usedforsecurity=False)
            ).hexdigest()
        if digest != metadata.md5sum:
            raise ValueError("ENCODE file does not match its published MD5")
        _encode_artifact(staged, requested)
        try:
            os.rename(staged, package_root)
        except OSError:
            if not package_root.is_dir():
                raise
    return _encode_artifact(package_root, requested)


ARTIFACT_PROVIDERS: dict[str, Callable[[IdentifierRef, str, Path], Artifact]] = {
    "ncbi_datasets": _refseq_gcf_path,
    "uniprot_rest": _uniprot_path,
    "ncbi_file": _file_path,
    "file_download": _file_path,
    "encode_file": _encode_path,
}


def artifact_capabilities() -> dict[str, Any]:
    """Describe supported artifacts without network or artifact-cache access.

    Returns:
        Independent JSON-compatible manifest containing its version and namespace
        records with provider names, defaults, and file representations.
    """
    return deepcopy(_artifact_capabilities())


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
        ValueError: If no provider supports the pair.
    """
    requested = (
        identifier
        if isinstance(identifier, IdentifierRef)
        else parse_identifier(identifier)
    )
    capability = _artifact_capabilities()["namespaces"].get(requested.namespace)
    if capability is None or (
        artifact is not None and artifact not in capability["kinds"]
    ):
        detail = f" and artifact {artifact!r}" if artifact is not None else ""
        raise ValueError(
            f"No artifact provider for namespace {requested.namespace!r}{detail}"
        )
    kind = capability["default_kind"] if artifact is None else artifact
    package_root = (
        settings.home
        / "artifacts"
        / _cache_component(requested.namespace, "namespace")
        / _cache_component(requested.accession, "accession")
    )
    return ARTIFACT_PROVIDERS[capability["provider"]](requested, kind, package_root)


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
    "artifact_capabilities",
    "open",
    "path",
]
