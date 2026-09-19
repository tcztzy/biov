"""Acceptance tests for the packaged identifiers.org registry response."""

import json
import re
from copy import deepcopy
from importlib.resources import files
from pathlib import Path

import pytest
from typer.testing import CliRunner

import biov.cli as cli
from biov.registry import (
    RegistryAssetError,
    namespaces_by_prefix,
    load_registry_asset,
    parse_registry_response,
    update_registry_asset,
)


def _packaged_response_bytes() -> bytes:
    """Read the packaged upstream response body."""
    return files("biov.assets").joinpath("identifiers_org_registry.json").read_bytes()


def test_v46_registry_asset_is_raw_upstream_response() -> None:
    """Keep the native response shape without a BioV storage protocol."""
    content = _packaged_response_bytes()
    response = json.loads(content)

    assert load_registry_asset() == response
    assert {
        "assetSchemaVersion",
        "source",
        "responseMetadata",
        "payloadMetadata",
        "resources",
        "institutions",
        "locations",
    }.isdisjoint(response)
    assert response["payload"]["namespaces"]
    assert all(
        isinstance(namespace["resources"], list)
        for namespace in response["payload"]["namespaces"]
    )


def test_every_namespace_has_a_compilable_runtime_rule() -> None:
    """Index only native namespace fields required by runtime routing."""
    namespaces = namespaces_by_prefix()
    source_namespaces = load_registry_asset()["payload"]["namespaces"]

    assert namespaces == {
        record["prefix"].casefold(): record for record in source_namespaces
    }
    assert all(re.compile(record["pattern"]) for record in namespaces.values())
    assert namespaces["refseq.gcf"]["id"] == 3719


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("id", "1"),
        ("prefix", 1),
        ("name", None),
        ("description", {}),
        ("pattern", "["),
        ("sampleId", 1),
        ("namespaceEmbeddedInLui", "false"),
        ("deprecated", 0),
    ],
)
def test_registry_rejects_invalid_runtime_fields(field: str, value: object) -> None:
    """Reject malformed upstream fields without silently coercing their types."""
    record = {**load_registry_asset()["payload"]["namespaces"][0], field: value}
    with pytest.raises(RegistryAssetError):
        parse_registry_response(json.dumps({"payload": {"namespaces": [record]}}))


def test_registry_asset_and_updater_are_packaged_project_assets() -> None:
    """Keep the raw snapshot importable and its updater version-controlled."""
    asset = files("biov.assets").joinpath("identifiers_org_registry.json")
    script = Path(__file__).parents[1] / "scripts" / "update_identifiers_registry.py"

    assert asset.is_file()
    assert script.is_file()


def test_registry_update_skips_a_byte_identical_response(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Leave a byte-identical raw asset untouched."""
    output = tmp_path / "registry.json"
    content = _packaged_response_bytes()
    output.write_bytes(content)
    inode = output.stat().st_ino

    monkeypatch.setattr(
        "biov.registry.fetch_registry",
        lambda **_kwargs: content,
    )
    result = update_registry_asset(output)

    assert result.status == "unchanged"
    assert output.stat().st_ino == inode
    assert output.read_bytes() == content


def test_registry_update_writes_upstream_bytes_without_restructuring(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Preserve unknown nested upstream data and the exact fetched body."""
    output = tmp_path / "registry.json"
    output.write_bytes(_packaged_response_bytes())
    response = deepcopy(load_registry_asset())
    response["futureMetadata"] = {"items": [1, {"native": True}]}
    response["payload"]["namespaces"][0]["futureNamespaceField"] = {"nested": ["value"]}
    content = json.dumps(response, ensure_ascii=False, separators=(", ", ": ")).encode()

    monkeypatch.setattr(
        "biov.registry.fetch_registry",
        lambda **_kwargs: content,
    )
    result = update_registry_asset(output)

    assert result.status == "updated"
    assert result.asset == response
    assert output.read_bytes() == content


def test_registry_update_rejects_invalid_raw_body_before_replacement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Validate the fetched body before atomically publishing it."""
    output = tmp_path / "registry.json"
    original = _packaged_response_bytes()
    output.write_bytes(original)
    invalid = b'{"payload":{"namespaces":"invalid"}}'

    monkeypatch.setattr(
        "biov.registry.fetch_registry",
        lambda **_kwargs: invalid,
    )

    with pytest.raises(RegistryAssetError, match="namespaces"):
        update_registry_asset(output)
    assert output.read_bytes() == original


def test_biov_registry_update_subcommand(monkeypatch, tmp_path: Path) -> None:
    """Expose raw registry synchronization through the BioV CLI."""
    output = tmp_path / "registry.json"
    expected_response = load_registry_asset()

    def fake_update(path, *, force, timeout):
        assert path == output
        assert force is False
        assert timeout == 60
        return cli.RegistryUpdateResult(
            status="unchanged",
            asset=expected_response,
            output=output,
        )

    monkeypatch.setattr(cli, "update_registry_asset", fake_update)
    result = CliRunner().invoke(
        cli.app,
        ["update-identifiers-registry", "--output", str(output)],
    )

    assert result.exit_code == 0
    assert "already current" in result.stdout


def test_biov_registry_update_failure_is_a_clean_cli_error(
    monkeypatch, tmp_path: Path
) -> None:
    """Map registry update failures to a stable exit status, not a traceback."""
    output = tmp_path / "registry.json"

    def failing_update(path, *, force, timeout):
        raise RegistryAssetError("registry response is not valid JSON")

    monkeypatch.setattr(cli, "update_registry_asset", failing_update)
    result = CliRunner().invoke(
        cli.app,
        ["update-identifiers-registry", "--output", str(output)],
    )

    assert result.exit_code == 2
    assert "not valid JSON" in result.stderr
    assert "Traceback" not in result.output
