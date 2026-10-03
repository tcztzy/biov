"""Validate the shared artifact capability asset."""

import json
from functools import cache
from importlib.resources import files
from typing import Any


def _unique_manifest_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Decode a manifest object without discarding duplicate keys.

    Returns:
        Decoded object with unique names.

    Raises:
        ValueError: If a namespace, kind, or other object key repeats.
    """
    value = dict(pairs)
    if len(value) != len(pairs):
        raise ValueError("artifact capabilities contain duplicate keys")
    return value


@cache
def load_capabilities() -> dict[str, Any]:
    """Load and validate the packaged artifact capability manifest.

    Returns:
        Manifest used for provider dispatch and package member selection.

    Raises:
        ValueError: If the manifest has an unsupported version or invalid record.
    """
    manifest = json.loads(
        files("biov.assets")
        .joinpath("artifact_capabilities.json")
        .read_text(encoding="utf-8"),
        object_pairs_hook=_unique_manifest_object,
    )
    if not isinstance(manifest, dict) or manifest.get("version") != 1:
        raise ValueError("artifact capabilities require version 1")
    namespaces = manifest.get("namespaces")
    if not isinstance(namespaces, dict) or not namespaces:
        raise ValueError("artifact capabilities require namespace records")
    for namespace, record in namespaces.items():
        if not namespace or not isinstance(record, dict):
            raise ValueError("artifact capabilities contain an invalid namespace")
        provider = record.get("provider")
        if not isinstance(provider, str) or provider not in {
            "ncbi_datasets",
            "uniprot_rest",
            "ncbi_file",
            "file_download",
            "encode_file",
        }:
            raise ValueError(f"unsupported artifact provider for {namespace!r}")
        kinds = record.get("kinds")
        default = record.get("default_kind")
        if (
            not isinstance(kinds, dict)
            or not kinds
            or not isinstance(default, str)
            or default not in kinds
        ):
            raise ValueError(f"invalid artifact kinds or default for {namespace!r}")
        for kind, properties in kinds.items():
            if not kind or not isinstance(properties, dict):
                raise ValueError(f"invalid artifact kind for {namespace!r}")
            if provider == "ncbi_datasets" and any(
                not isinstance(properties.get(field), str) or not properties[field]
                for field in ("fileType", "label")
            ):
                raise ValueError("NCBI artifact kinds require fileType and label")
    return manifest
