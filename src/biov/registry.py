"""Load, validate, index, and refresh the native identifiers.org registry."""

import json
import os
import re
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Literal, TypedDict, cast
from urllib.parse import quote
from urllib.request import Request, urlopen

from pydantic import TypeAdapter, ValidationError

REGISTRY_DATASET_URL = (
    "https://registry.api.identifiers.org/resolutionApi/getResolverDataset"
)
REGISTRY_ASSET_NAME = "identifiers_org_registry.json"
DATA_RESOURCE_NAMESPACE_PREFIXES = frozenset({"refseq.gcf", "uniprot"})
URI_ACCESSION_SAFE = "/:;,@!$&'*+=~"

RegistryResponse = dict[str, Any]


class RegistryAssetError(ValueError):
    """The registry response lacks data required by the runtime."""


class Namespace(TypedDict):
    """Required fields of an unmodified upstream namespace record."""

    id: int
    prefix: str
    name: str
    description: str
    pattern: str
    sampleId: str | None
    namespaceEmbeddedInLui: bool
    deprecated: bool


@dataclass(frozen=True, slots=True)
class RegistryUpdateResult:
    """Outcome of synchronizing a local registry response."""

    status: Literal["updated", "unchanged"]
    asset: RegistryResponse
    output: Path


_NAMESPACE_SCHEMA = TypeAdapter(list[Namespace])


def parse_registry_response(content: bytes | str) -> RegistryResponse:
    """Validate runtime fields while preserving the complete native response.

    Returns:
        Parsed upstream JSON, including unknown fields and nesting.

    Raises:
        RegistryAssetError: If JSON or required runtime fields are invalid.
    """
    try:
        response = json.loads(content)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RegistryAssetError("registry response is not valid JSON") from error
    if not isinstance(response, dict):
        raise RegistryAssetError("registry response must be a JSON object")
    try:
        records = _NAMESPACE_SCHEMA.validate_python(
            response["payload"]["namespaces"], strict=True
        )
    except (KeyError, TypeError, ValidationError) as error:
        raise RegistryAssetError(
            f"invalid registry response.payload.namespaces: {error}"
        ) from error
    seen: set[str] = set()
    for record in records:
        prefix = record["prefix"].casefold()
        if not prefix or prefix in seen:
            raise RegistryAssetError("namespace prefixes must be non-empty and unique")
        try:
            re.compile(record["pattern"])
        except re.error as error:
            raise RegistryAssetError(
                f"invalid pattern for namespace {prefix!r}"
            ) from error
        seen.add(prefix)
    return response


def load_registry_asset() -> RegistryResponse:
    """Read and validate the complete packaged upstream response.

    Returns:
        The complete native JSON response.
    """
    return parse_registry_response(
        files("biov.assets").joinpath(REGISTRY_ASSET_NAME).read_bytes()
    )


@cache
def namespaces_by_prefix() -> dict[str, Namespace]:
    """Index native records once for both identifier parsing and MCP resources.

    Returns:
        Native records keyed by case-insensitive prefix.
    """
    records = cast("list[Namespace]", load_registry_asset()["payload"]["namespaces"])
    return {record["prefix"].casefold(): record for record in records}


def build_namespace_resource_uri(namespace: Namespace, accession: str) -> str:
    """Build a percent-safe data or generic identifiers resource URI.

    Returns:
        The data resource URI or generic identifiers resource URI.
    """
    encoded = quote(accession, safe=URI_ACCESSION_SAFE)
    prefix = namespace["prefix"]
    if prefix in DATA_RESOURCE_NAMESPACE_PREFIXES:
        return f"{prefix}://{encoded}"
    return f"identifiers://{quote(prefix, safe='._~-')}:{encoded}"


def default_registry_asset_path() -> Path:
    """Return the writable asset path in an unpacked BioV installation.

    Raises:
        RegistryAssetError: If the package is not writable in place.
    """
    asset = files("biov.assets").joinpath(REGISTRY_ASSET_NAME)
    if not isinstance(asset, Path):
        raise RegistryAssetError(
            "the packaged registry is not a writable filesystem asset; pass --output"
        )
    return asset


def fetch_registry(*, timeout: float = 60) -> bytes:
    """Fetch the complete official resolver response without transforming it.

    Args:
        timeout: Network timeout in seconds.

    Returns:
        Unmodified HTTP response body.

    Raises:
        ValueError: If the timeout is not positive.
    """
    if timeout <= 0:
        raise ValueError("timeout must be greater than zero")
    request = Request(
        REGISTRY_DATASET_URL,
        headers={
            "Accept": "application/json",
            "User-Agent": "BioV identifiers.org registry updater",
        },
    )
    with urlopen(request, timeout=timeout) as response:  # noqa: S310
        return response.read()


def _write_asset_atomically(content: bytes, output: Path) -> None:
    """Write validated upstream bytes without exposing a partial file."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent, prefix=f".{output.name}.") as staging:
        temporary_path = Path(staging) / output.name
        with temporary_path.open("wb") as temporary:
            temporary.write(content)
            temporary.flush()
            os.fsync(temporary.fileno())
        os.chmod(temporary_path, 0o644)
        os.replace(temporary_path, output)
    namespaces_by_prefix.cache_clear()


def update_registry_asset(
    output: Path | None = None,
    *,
    force: bool = False,
    timeout: float = 60,
) -> RegistryUpdateResult:
    """Replace a registry asset with the unmodified upstream response body.

    Args:
        output: Destination JSON file; defaults to the packaged asset.
        force: Replace the destination even when its bytes are unchanged.
        timeout: Network timeout in seconds.

    Returns:
        Update status, native parsed response, and destination path.
    """
    destination = output or default_registry_asset_path()
    try:
        existing = destination.read_bytes()
    except FileNotFoundError:
        existing = None
    content = fetch_registry(timeout=timeout)
    response = parse_registry_response(content)
    if not force and existing == content:
        return RegistryUpdateResult("unchanged", response, destination)
    _write_asset_atomically(content, destination)
    return RegistryUpdateResult("updated", response, destination)


__all__ = [
    "DATA_RESOURCE_NAMESPACE_PREFIXES",
    "REGISTRY_ASSET_NAME",
    "REGISTRY_DATASET_URL",
    "URI_ACCESSION_SAFE",
    "Namespace",
    "RegistryAssetError",
    "RegistryResponse",
    "RegistryUpdateResult",
    "build_namespace_resource_uri",
    "default_registry_asset_path",
    "fetch_registry",
    "load_registry_asset",
    "namespaces_by_prefix",
    "parse_registry_response",
    "update_registry_asset",
]
