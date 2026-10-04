"""Acceptance tests for the identifiers.org MCP interface."""

import io
import json
import tomllib
from collections.abc import Awaitable, Callable
from pathlib import Path
from typing import Any, cast

import anyio
import pytest
from mcp.server.lowlevel.helper_types import ReadResourceContents
from mcp.server.mcpserver.exceptions import (
    ResourceError,
    ToolError,
    UnexpectedResourceError,
)
from typer.testing import CliRunner

import biov.artifacts as artifacts_module
import biov.identifiers as identifiers_module
import biov.mcp as mcp_module
from biov import artifact_capabilities, cli
from biov.artifacts import Artifact
from biov.identifiers import (
    MAX_PROMPT_CANDIDATES,
    IdentifierNotFoundError,
    IdentifierServiceError,
    extract_identifier_candidates,
    parse_identifier,
)
from biov.mcp import create_mcp_server
from biov.registry import namespaces_by_prefix

Resolution = dict[str, Any]
Resolver = Callable[[str], Awaitable[Resolution]]


def _resolution(
    compact_id: str,
    namespace: str,
    local_id: str,
    provider_code: str | None = None,
) -> Resolution:
    """Build the resolver fragment used by prompt-link tests.

    Returns:
        Minimal resolver response naming the requested namespace and local ID.
    """
    return {
        "payload": {
            "parsedCompactIdentifier": {
                "providerCode": provider_code,
                "namespace": namespace,
                "localId": local_id,
                "rawRequest": compact_id,
            },
            "resolvedResources": [],
        }
    }


def test_resolver_returns_native_json_and_rejects_malformed_responses(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep the default async resolver's HTTP boundary intact."""
    expected = _resolution("uniprot:P12345", "uniprot", "P12345")
    body = json.dumps(expected).encode()

    def urlopen(request, *, timeout):
        assert request.full_url == "https://resolver.api.identifiers.org/uniprot:P12345"
        assert timeout == 10
        return io.BytesIO(body)

    monkeypatch.setattr(identifiers_module, "urlopen", urlopen)
    assert (
        anyio.run(identifiers_module.resolve_identifier, "uniprot:P12345") == expected
    )
    body = b"not JSON"
    with pytest.raises(IdentifierServiceError, match="invalid JSON"):
        anyio.run(identifiers_module.resolve_identifier, "uniprot:P12345")


def test_extract_identifier_candidates_preserves_order_and_scope() -> None:
    """Parse explicit identifiers without guessing bare accessions."""
    prompt = (
        "Compare uniprot:P12345, "
        "https://identifiers.org/ols/TAXONOMY:9606; "
        "GO:0006915 and doi:10.1038/s41586-020-2649-2. "
        "Ignore bare P12345 and repeat uniprot:P12345."
    )

    assert extract_identifier_candidates(prompt) == [
        "uniprot:P12345",
        "ols/TAXONOMY:9606",
        "GO:0006915",
        "doi:10.1038/s41586-020-2649-2",
    ]


def test_extract_identifier_candidates_accepts_identifiers_resource_uri() -> None:
    """Map the generic identifiers scheme back to resolver input."""
    assert extract_identifier_candidates(
        "Read identifiers://doi:10.1038/s41586-020-2649-2",
    ) == ["doi:10.1038/s41586-020-2649-2"]


@pytest.mark.parametrize(
    "reference",
    [
        "refseq.gcf://GCF_000001030.2",
        "GCF_000001030.2",
        "UNIPROT://P42212",
        "identifiers://go:GO%3A0006915",
        "identifiers://doi:10.1038%2Fs41586-020-2649-2",
        "identifiers://3dmet:B00162",
    ],
)
def test_prompt_and_single_identifier_share_resource_and_bare_rules(
    reference: str,
) -> None:
    """Keep prompt variants aligned with the single-identifier parser."""
    assert extract_identifier_candidates(f"Read `{reference}`.") == [
        parse_identifier(reference).compact_id
    ]


def test_extract_identifier_candidates_is_bounded() -> None:
    """Bound work induced by an arbitrarily long prompt."""
    prompt = " ".join(f"example{i}:value" for i in range(MAX_PROMPT_CANDIDATES + 5))

    assert len(extract_identifier_candidates(prompt)) == MAX_PROMPT_CANDIDATES


def test_parse_identifiers_uses_registry_uris_and_deduplicates_links() -> None:
    """Map resolver-parsed namespaces to generated, location-independent links."""
    calls: list[str] = []
    resolutions = {
        "refseq.gcf:GCF_000001030.2": _resolution(
            "refseq.gcf:GCF_000001030.2",
            "refseq.gcf",
            "GCF_000001030.2",
        ),
        "ols/taxonomy:9606": _resolution(
            "ols/taxonomy:9606",
            "taxonomy",
            "9606",
            "ols",
        ),
        "taxonomy:9606": _resolution("taxonomy:9606", "taxonomy", "9606"),
        "GO:0006915": _resolution("GO:0006915", "go", "GO:0006915"),
    }

    async def resolve(compact_id: str) -> Resolution:
        calls.append(compact_id)
        if compact_id == "invalid:1":
            raise IdentifierNotFoundError(compact_id, "unknown namespace")
        return resolutions[compact_id]

    server = create_mcp_server(resolve)

    async def call_tool() -> Any:
        return await server.call_tool(
            "parse_identifiers",
            {
                "prompt": (
                    "refseq.gcf:GCF_000001030.2 invalid:1 "
                    "ols/taxonomy:9606 taxonomy:9606 GO:0006915"
                )
            },
        )

    result = anyio.run(call_tool)
    links = [content for content in result.content if content.type == "resource_link"]

    assert calls == [
        "refseq.gcf:GCF_000001030.2",
        "invalid:1",
        "ols/taxonomy:9606",
        "taxonomy:9606",
        "GO:0006915",
    ]
    assert [str(link.uri) for link in links] == [
        "refseq.gcf://GCF_000001030.2",
        "identifiers://taxonomy:9606",
        "identifiers://go:GO:0006915",
    ]
    assert json.loads(result.content[0].text) == {
        "count": 3,
        "resource_uris": [str(link.uri) for link in links],
    }


@pytest.mark.parametrize(
    "reference",
    [
        "refseq.gcf://GCF_000001030.2",
        "identifiers://refseq.gcf:GCF_000001030.2",
        "refseq.gcf:GCF_000001030.2",
        "GCF_000001030.2",
    ],
)
def test_parse_identifiers_accepts_each_gcf_variant_in_isolation(
    reference: str,
) -> None:
    """Recognize each common GCF form without relying on another valid form."""
    calls: list[str] = []

    async def resolve(compact_id: str) -> Resolution:
        calls.append(compact_id)
        return _resolution(
            compact_id,
            "refseq.gcf",
            "GCF_000001030.2",
        )

    server = create_mcp_server(resolve)

    async def call_tool() -> Any:
        return await server.call_tool(
            "parse_identifiers",
            {
                "prompt": (
                    f"请问 `{reference}` 的 GC content 是多少？ "
                    "Do not infer GCA_000155495.1 or P12345."
                )
            },
        )

    result = anyio.run(call_tool)
    links = [content for content in result.content if content.type == "resource_link"]

    assert calls == ["refseq.gcf:GCF_000001030.2"]
    assert [str(link.uri) for link in links] == ["refseq.gcf://GCF_000001030.2"]


def test_parse_identifiers_deduplicates_gcf_variants_before_resolver() -> None:
    """Canonicalize equivalent URI, Compact, and bare forms before resolution."""
    calls: list[str] = []

    async def resolve(compact_id: str) -> Resolution:
        calls.append(compact_id)
        return _resolution(compact_id, "refseq.gcf", "GCF_000001030.2")

    server = create_mcp_server(resolve)

    async def call_tool() -> Any:
        return await server.call_tool(
            "parse_identifiers",
            {
                "prompt": (
                    "refseq.gcf://GCF_000001030.2, "
                    "refseq.gcf:GCF_000001030.2 and GCF_000001030.2"
                )
            },
        )

    result = anyio.run(call_tool)
    links = [content for content in result.content if content.type == "resource_link"]

    assert calls == ["refseq.gcf:GCF_000001030.2"]
    assert len(links) == 1


def test_parse_identifiers_rejects_invalid_bare_gcf_before_resolver() -> None:
    """Apply the allowlisted registry rule before resolving shorthand input."""
    calls: list[str] = []

    async def resolve(compact_id: str) -> Resolution:
        calls.append(compact_id)
        return _resolution(compact_id, "refseq.gcf", compact_id)

    server = create_mcp_server(resolve)

    async def call_tool() -> Any:
        return await server.call_tool(
            "parse_identifiers",
            {"prompt": "What is the GC content of GCF_00001030?"},
        )

    result = anyio.run(call_tool)

    assert calls == []
    assert json.loads(result.content[0].text) == {"count": 0, "resource_uris": []}


def test_parse_identifiers_surfaces_service_failure() -> None:
    """Do not misreport an upstream failure as an absent identifier."""

    async def unavailable(compact_id: str) -> Resolution:
        raise IdentifierServiceError("identifiers.org service unavailable")

    server = create_mcp_server(unavailable)

    async def call_tool() -> Any:
        return await server.call_tool(
            "parse_identifiers",
            {"prompt": "uniprot:P12345"},
        )

    with pytest.raises(ToolError, match=r"identifiers\.org service unavailable"):
        anyio.run(call_tool)


def test_identifier_resource_returns_native_resolution_with_reserved_id() -> None:
    """Return resolver JSON while retaining slash-bearing local IDs."""
    seen: list[str] = []

    async def resolve(compact_id: str) -> Resolution:
        seen.append(compact_id)
        return {"compactIdentifier": compact_id}

    server = create_mcp_server(resolve)

    async def read_resource() -> Any:
        return await server.read_resource("identifiers://doi:10.1038/s41586-020-2649-2")

    contents = list(anyio.run(read_resource))

    assert seen == ["doi:10.1038/s41586-020-2649-2"]
    assert contents[0].mime_type == "application/json"
    assert json.loads(contents[0].content) == {
        "compactIdentifier": "doi:10.1038/s41586-020-2649-2"
    }


def test_identifier_resource_handles_embedded_prefix() -> None:
    """Do not duplicate namespace prefixes embedded in local identifiers."""
    seen: list[str] = []

    async def resolve(compact_id: str) -> Resolution:
        seen.append(compact_id)
        return {"compactIdentifier": compact_id}

    server = create_mcp_server(resolve)

    async def read_resource() -> Any:
        return await server.read_resource("identifiers://go:GO:0006915")

    anyio.run(read_resource)

    assert seen == ["GO:0006915"]


def test_identifier_resource_supports_numeric_registry_prefix() -> None:
    """Use the generic scheme for prefixes that cannot be URI schemes."""
    seen: list[str] = []

    async def resolve(compact_id: str) -> Resolution:
        seen.append(compact_id)
        return {"compactIdentifier": compact_id}

    server = create_mcp_server(resolve)

    async def read_resource() -> Any:
        return await server.read_resource("identifiers://3dmet:B00162")

    anyio.run(read_resource)

    assert seen == ["3dmet:B00162"]


def test_namespace_resource_rejects_accession_outside_registry_pattern() -> None:
    """Apply the packaged namespace rule before making a resolver request."""
    seen: list[str] = []

    async def resolve(compact_id: str) -> Resolution:
        seen.append(compact_id)
        return {"compactIdentifier": compact_id}

    server = create_mcp_server(resolve)

    async def read_resource() -> Any:
        return await server.read_resource("refseq.gcf://GCF_00001030")

    with pytest.raises(ResourceError, match=r"does not match registry pattern"):
        anyio.run(read_resource)
    assert seen == []


def test_registry_resource_returns_native_namespace_record() -> None:
    """Expose one losslessly reconstructed identifiers.org namespace object."""

    async def unexpected_resolver(compact_id: str) -> Resolution:
        raise AssertionError(f"registry read contacted resolver for {compact_id}")

    server = create_mcp_server(unexpected_resolver)

    async def read_resource() -> Any:
        return await server.read_resource("identifiers://go")

    contents = list(anyio.run(read_resource))

    assert contents[0].mime_type == "application/json"
    assert json.loads(contents[0].content) == namespaces_by_prefix()["go"]


@pytest.mark.parametrize(
    ("namespace", "accession"),
    [
        ("pubmed", "16333295"),
        ("clinvar", "12345"),
        ("dbsnp", "rs121909098"),
        ("geo", "GDS1234"),
    ],
)
def test_ncbi_identifier_forms_produce_direct_resource_links(
    namespace: str, accession: str
) -> None:
    """Route each explicit ID form through the same data resource."""
    compact_id = f"{namespace}:{accession}"
    uri = f"{namespace}://{accession}"
    seen: list[str] = []

    async def resolve(value: str) -> Resolution:
        seen.append(value)
        return _resolution(value, namespace, accession)

    server = create_mcp_server(resolve)
    for value in (compact_id, uri, f"https://identifiers.org/{compact_id}"):
        assert parse_identifier(value).resource_uri == uri

        async def call_tool(value: str = value) -> Any:
            return await server.call_tool("parse_identifiers", {"prompt": value})

        result = anyio.run(call_tool)
        assert str(result.content[1].uri) == uri
        assert result.content[1].meta["resourceKind"] == "data"
    assert seen == [compact_id] * 3


@pytest.mark.parametrize(
    "uri", ["pubmed://PMID123", "clinvar://VCV123", "dbsnp://123", "geo://GEX123"]
)
def test_ncbi_resources_validate_registry_rule_before_io(
    uri: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reject invalid local IDs before contacting NCBI."""

    def unexpected_request(*args, **kwargs):
        raise AssertionError("Invalid accession reached NCBI")

    monkeypatch.setattr(mcp_module, "path", unexpected_request)
    server = create_mcp_server()

    async def read_resource() -> Any:
        return await server.read_resource(uri)

    with pytest.raises(ResourceError, match="does not match registry pattern"):
        anyio.run(read_resource)


def _artifact_path(
    tmp_path: Path,
    calls: list[tuple[str, str]],
) -> Callable[..., Artifact]:
    """Create deterministic raw provider files for MCP data-resource tests.

    Returns:
        Callable writing the requested raw file and recording the request.
    """

    def path(identifier: str, *, artifact: str) -> Artifact:
        reference = parse_identifier(identifier)
        calls.append((identifier, artifact))
        package_root = tmp_path / reference.namespace / reference.accession
        package_root.mkdir(parents=True, exist_ok=True)
        path = package_root / f"{reference.accession}.json"
        path.write_text('{\n  "primaryAccession": "P42212"\n}\n')
        return Artifact(
            path=path,
            package_root=package_root,
            requested_identifier=reference,
            identifier=reference,
            kind=artifact,
            size=path.stat().st_size,
        )

    return path


def test_uniprot_resource_returns_original_provider_json(
    tmp_path: Path,
) -> None:
    """Return cached provider-native JSON text without a BioV envelope."""
    calls: list[tuple[str, str]] = []
    server = create_mcp_server(artifact_path=_artifact_path(tmp_path, calls))

    async def read_resource() -> Any:
        return await server.read_resource("uniprot://P42212")

    contents = list(anyio.run(read_resource))

    assert contents[0].mime_type == "application/json"
    assert contents[0].content == '{\n  "primaryAccession": "P42212"\n}\n'
    assert calls == [("uniprot://P42212", "entry_json")]


def test_v47_refseq_resource_returns_raw_summary_without_artifact_path() -> None:
    """Keep MCP metadata reads independent from complete analysis packages."""
    summary = (
        '{"reports":[{"accession":"GCF_000001030.2",'
        '"organism":{"organism_name":"Brucella abortus"}}],"total_count":1}\n'
    )
    calls: list[str] = []

    def read_summary(accession: str) -> str:
        calls.append(accession)
        return summary

    def unexpected_artifact_path(identifier: str, *, artifact: str) -> Artifact:
        raise AssertionError(f"unexpected path resolution: {identifier} {artifact}")

    server = create_mcp_server(
        artifact_path=unexpected_artifact_path,
        refseq_reader=read_summary,
    )

    async def read_resource() -> Any:
        return await server.read_resource("refseq.gcf://GCF_000001030")

    contents = list(anyio.run(read_resource))

    assert contents[0].mime_type == "application/json"
    assert contents[0].content == summary
    assert calls == ["GCF_000001030"]


@pytest.mark.parametrize(
    ("uri", "failure"),
    [
        (
            "uniprot://P42212",
            OSError("Incomplete HTTP response: received 0 of 971 bytes"),
        ),
        (
            "refseq.gcf://GCF_000001030",
            ValueError("NCBI genome summary must be a JSON object"),
        ),
    ],
)
def test_provider_failures_reach_the_client_with_their_message(
    uri: str, failure: Exception
) -> None:
    """Classify every provider failure so the SDK does not withhold its text."""

    def fail_path(identifier: str, *, artifact: str) -> Artifact:
        raise failure

    def fail_reader(accession: str) -> str:
        raise failure

    server = create_mcp_server(artifact_path=fail_path, refseq_reader=fail_reader)

    async def read_resource() -> Any:
        return await server.read_resource(uri)

    with pytest.raises(ResourceError, match=str(failure)) as raised:
        anyio.run(read_resource)
    assert not isinstance(raised.value, UnexpectedResourceError)


def test_uniprot_resource_fetches_and_caches_only_raw_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reuse the REST JSON cache without requesting the optional FASTA."""
    downloads: list[str] = []

    def download_json(accession: str, destination: Path) -> None:
        """Reject any request for the optional FASTA representation."""
        assert destination.suffix == ".json"
        downloads.append(accession)
        destination.write_text('{\n  "primaryAccession": "P42212"\n}\n')

    monkeypatch.setattr(artifacts_module, "download_uniprot_entry", download_json)
    monkeypatch.setattr(artifacts_module.settings, "home", tmp_path)
    server = create_mcp_server()

    async def read_twice() -> list[str]:
        first = cast(
            "ReadResourceContents",
            next(iter(await server.read_resource("uniprot://P42212"))),
        )
        second = cast(
            "ReadResourceContents",
            next(iter(await server.read_resource("uniprot://P42212"))),
        )
        assert isinstance(first.content, str)
        assert isinstance(second.content, str)
        return [first.content, second.content]

    assert anyio.run(read_twice) == [
        '{\n  "primaryAccession": "P42212"\n}\n',
        '{\n  "primaryAccession": "P42212"\n}\n',
    ]
    assert downloads == ["P42212"]
    package_root = tmp_path / "artifacts" / "uniprot" / "P42212"
    assert (package_root / "P42212.json").is_file()
    assert not (package_root / "P42212.fasta").exists()


def test_resolve_identifiers_tool_embeds_the_same_resource() -> None:
    """Support clients that can call tools but cannot issue resources/read."""

    async def resolve(compact_id: str) -> Resolution:
        return {"compactIdentifier": compact_id}

    server = create_mcp_server(resolve)
    uri = "identifiers://doi:10.1038/s41586-020-2649-2"

    async def call_tool() -> Any:
        return await server.call_tool("resolve_identifiers", {"uri": uri})

    result = anyio.run(call_tool)
    embedded = result.content[0]

    assert embedded.type == "resource"
    assert str(embedded.resource.uri) == uri
    assert embedded.resource.mime_type == "application/json"
    assert json.loads(embedded.resource.text) == {
        "compactIdentifier": "doi:10.1038/s41586-020-2649-2"
    }


def test_mcp_surface_has_file_and_identifier_templates() -> None:
    """Publish supported files and identifier resources without catalogs."""

    async def resolve(compact_id: str) -> Resolution:
        return {"compactIdentifier": compact_id}

    server = create_mcp_server(resolve)
    templates = anyio.run(server.list_resource_templates)
    tools = anyio.run(server.list_tools)
    by_uri = {template.uri_template: template for template in templates}

    assert set(by_uri) == {
        *(
            f"{namespace}://{{+accession}}"
            for namespace in artifact_capabilities()["namespaces"]
        ),
        "identifiers://{registry}",
        "identifiers://{registry}:{+id}",
        "file://{+path}",
    }
    refseq = by_uri["refseq.gcf://{+accession}"]
    assert refseq.meta == {
        "namespaceId": 3719,
        "namespacePrefix": "refseq.gcf",
        "identifierPattern": r"^GCF_[0-9]{9}(\.[0-9]+)?$",
        "sampleId": "GCF_000001405",
    }
    assert {tool.name for tool in tools} == {
        "parse_identifiers",
        "resolve_identifiers",
        "run_analysis",
        "inspect_analysis",
    }
    for tool in (
        t for t in tools if t.name in {"parse_identifiers", "resolve_identifiers"}
    ):
        annotations = tool.annotations
        assert annotations is not None
        assert annotations.read_only_hint is True
        assert annotations.destructive_hint is False
        assert annotations.idempotent_hint is True
        assert annotations.open_world_hint is True


def test_biov_mcp_handler_is_preserved_with_native_entry_point(
    monkeypatch,
) -> None:
    """Launch MCP through ``biov mcp`` without a standalone console script."""
    started: list[bool] = []
    monkeypatch.setattr("biov.mcp.main", lambda: started.append(True))

    result = CliRunner().invoke(cli.app, ["mcp"])

    assert result.exit_code == 0
    assert started == [True]

    pyproject_path = Path(__file__).parents[1] / "pyproject.toml"
    with pyproject_path.open("rb") as file:
        pyproject = tomllib.load(file)

    scripts = pyproject["project"]["scripts"]
    assert "biov" not in scripts
    assert [
        binary["target"] for binary in pyproject["tool"]["setuptools-rust"]["bins"]
    ] == ["biov"]
    assert "biov-mcp" not in scripts


@pytest.mark.parametrize(
    ("uri", "kind"),
    [
        ("pubmed://23193287", "article_pdf"),
        ("clinvar://9", "variant_xml"),
        ("dbsnp://rs328", "refsnp_json"),
        ("geo://GSE100", "expression_matrix"),
        ("pdb://1CRN", "structure_cif"),
    ],
)
def test_file_resource_describes_available_data_without_downloading(
    uri, kind, monkeypatch
):
    """Describe provider representations without downloading any data."""

    def fail(*args, **kwargs):
        raise AssertionError("MCP file descriptions must not download data")

    server = create_mcp_server(artifact_path=fail)

    async def read():
        return list(await server.read_resource(uri))

    content = anyio.run(read)[0]
    assert isinstance(content, ReadResourceContents)
    result = json.loads(content.content)
    assert result["uri"] == uri
    assert result["default_kind"] == kind
    assert kind in result["representations"]
    assert "url" not in result["representations"][kind]


def test_provider_errors_reach_mcp_client(tmp_path):
    """Check wire responses, not the server process's private exception cause."""
    import sys

    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client
    from mcp.shared.exceptions import MCPError

    script = tmp_path / "failing_server.py"
    script.write_text("""from biov.mcp import create_mcp_server
async def unavailable(identifier):
    raise OSError("provider unavailable")
def missing(*args, **kwargs):
    raise FileNotFoundError("provider entry missing")
create_mcp_server(resolver=unavailable, artifact_path=missing).run()
""")

    async def check():
        with anyio.fail_after(30):
            async with stdio_client(
                StdioServerParameters(command=sys.executable, args=[str(script)])
            ) as (reader, writer):
                async with ClientSession(reader, writer) as session:
                    await session.initialize()
                    with pytest.raises(MCPError, match="provider entry missing"):
                        await session.read_resource("uniprot://P12345")
                    for name, arguments, message in (
                        (
                            "parse_identifiers",
                            {"prompt": "uniprot:P12345"},
                            "provider unavailable",
                        ),
                        (
                            "resolve_identifiers",
                            {"uri": "uniprot://P12345"},
                            "provider entry missing",
                        ),
                    ):
                        result = await session.call_tool(name, arguments)
                        assert result.is_error
                        assert result.content[0].type == "text"
                        assert message in result.content[0].text

    anyio.run(check)
