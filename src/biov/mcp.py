"""BioV Model Context Protocol server."""

import json
import re
from collections.abc import Awaitable, Callable
from typing import cast

from anyio.to_thread import run_sync
from mcp.server import MCPServer
from mcp.server.lowlevel.helper_types import ReadResourceContents
from mcp.server.mcpserver.exceptions import ResourceError, ToolError
from mcp.types import (
    EmbeddedResource,
    ResourceLink,
    TextContent,
    TextResourceContents,
    ToolAnnotations,
)

from .artifacts import (
    Artifact,
    ArtifactError,
    genome_summary,
    path,
)
from .identifiers import (
    BARE_IDENTIFIER_NAMESPACE_ALLOWLIST,
    IdentifierNotFoundError,
    IdentifierServiceError,
    Resolution,
    build_identifiers_org_url,
    extract_identifier_candidates,
    resolve_identifier,
)
from .registry import (
    DATA_RESOURCE_NAMESPACE_PREFIXES,
    Namespace,
    build_namespace_resource_uri,
    namespaces_by_prefix,
)

Resolver = Callable[[str], Awaitable[Resolution]]
ArtifactPath = Callable[..., Artifact]
DataReader = Callable[[str], str]


def _json_text(value: object) -> str:
    """Serialize one native upstream JSON value for MCP text transport.

    Returns:
        Indented UTF-8 JSON text without changing its data model.
    """
    return json.dumps(value, ensure_ascii=False, indent=2)


def _compact_identifier(
    namespace: Namespace,
    accession: str,
) -> str:
    """Reconstruct resolver input for one namespace resource.

    Returns:
        Resolver-compatible Compact Identifier.
    """
    if namespace["namespaceEmbeddedInLui"]:
        return accession
    return f"{namespace['prefix']}:{accession}"


def _validate_accession(namespace: Namespace, accession: str) -> None:
    """Apply one packaged identifiers.org accession rule.

    Raises:
        ResourceError: If the accession does not match the registry rule.
    """
    if re.fullmatch(namespace["pattern"], accession) is None:
        raise ResourceError(
            f"Accession {accession!r} does not match registry pattern "
            f"{namespace['pattern']!r} for {namespace['prefix']!r}"
        )


def _data_resource_handler(
    namespace: Namespace,
    reader: DataReader,
) -> Callable[[str], Awaitable[str]]:
    """Create one canonical provider-data resource handler.

    Returns:
        Asynchronous handler for one provider accession.
    """

    async def read_data(accession: str) -> str:
        """Return the provider-native canonical JSON representation.

        Raises:
            ResourceError: If validation or artifact access fails.
        """
        _validate_accession(namespace, accession)
        try:
            return await run_sync(reader, accession)
        except ArtifactError as error:
            raise ResourceError(str(error)) from error

    return read_data


def _parsed_identifier(resolution: Resolution) -> tuple[str, str]:
    """Read namespace and local ID from an official resolver response.

    Returns:
        Parsed namespace prefix and namespace-local identifier.

    Raises:
        IdentifierServiceError: If parsed identifier fields are absent or invalid.
    """
    payload = resolution.get("payload")
    parsed = (
        payload.get("parsedCompactIdentifier") if isinstance(payload, dict) else None
    )
    if not isinstance(parsed, dict):
        raise IdentifierServiceError(
            "identifiers.org resolver response has no parsed Compact Identifier"
        )
    namespace = parsed.get("namespace")
    local_id = parsed.get("localId")
    if not isinstance(namespace, str) or not isinstance(local_id, str):
        raise IdentifierServiceError(
            "identifiers.org resolver returned invalid parsed identifier fields"
        )
    return namespace, local_id


def create_mcp_server(
    resolver: Resolver | None = None,
    artifact_path: ArtifactPath | None = None,
    refseq_reader: DataReader | None = None,
) -> MCPServer:
    """Create the BioV MCP server.

    Args:
        resolver: Optional resolver implementation, primarily for tests and embedding.
        artifact_path: Optional artifact path function for data-backed resources.
        refseq_reader: Optional NCBI genome-summary reader.

    Returns:
        Configured server with data and identifiers.org resources plus tools.
    """
    resolve = resolver or resolve_identifier
    resource_path = artifact_path or path
    read_refseq = refseq_reader or genome_summary
    namespaces = namespaces_by_prefix()
    server = MCPServer(
        name="biov",
        title="BioV",
        description="LLM-native molecular biology resources",
        instructions=(
            "Call parse_identifiers with prompts containing biological IDs. "
            "Read refseq.gcf:// and uniprot:// resources for canonical upstream "
            "data. Read identifiers://<registry> for native registry metadata and "
            "identifiers://<registry>:<id> for native resolver JSON. Tool-only "
            "clients can pass any returned resource URI to resolve_identifiers. "
            "Generated Python should pass biov.path(uri) to libraries accepting "
            "os.PathLike, or use biov.open(uri, mode='rb' or 'rt') when a readable "
            "file object is required. Call both only inside the execution environment."
        ),
    )

    def read_uniprot_json(accession: str) -> str:
        """Read the original cached UniProt entry JSON.

        Returns:
            Original provider JSON text.

        Raises:
            ArtifactError: If reading the artifact fails.
        """
        artifact = resource_path(f"uniprot://{accession}", artifact="entry_json")
        try:
            return artifact.path.read_text(encoding="utf-8")
        except (OSError, UnicodeError) as error:
            raise ArtifactError("Could not read uniprot artifact JSON") from error

    for prefix, description, reader in (
        (
            "refseq.gcf",
            (
                "Return the original NCBI Datasets genome-summary JSON for a RefSeq "
                "GCF assembly without downloading its data package."
            ),
            read_refseq,
        ),
        (
            "uniprot",
            "Return the complete original UniProtKB REST JSON for one accession.",
            read_uniprot_json,
        ),
    ):
        namespace = namespaces[prefix]
        server.resource(
            f"{namespace['prefix']}://{{+accession}}",
            name=f"read_{namespace['prefix']}",
            title=namespace["name"],
            description=description,
            mime_type="application/json",
            meta={
                "namespaceId": namespace["id"],
                "namespacePrefix": namespace["prefix"],
                "identifierPattern": namespace["pattern"],
                "sampleId": namespace["sampleId"],
            },
        )(_data_resource_handler(namespace, reader))

    @server.resource(
        "identifiers://{registry}:{+id}",
        name="resolve_identifier_resource",
        title="Identifiers.org resolver response",
        description=(
            "Return the complete native identifiers.org resolver JSON for one "
            "registry-local identifier."
        ),
        mime_type="application/json",
    )
    async def read_identifier(registry: str, id: str) -> str:
        """Return the official resolver JSON for one registry-local ID.

        Raises:
            ResourceError: If validation or identifiers.org resolution fails.
        """
        namespace = namespaces.get(registry.casefold())
        if namespace is None:
            raise ResourceError(f"Unknown identifiers.org registry {registry!r}")
        _validate_accession(namespace, id)
        try:
            resolution = await resolve(_compact_identifier(namespace, id))
        except (IdentifierNotFoundError, IdentifierServiceError) as error:
            raise ResourceError(str(error)) from error
        return _json_text(resolution)

    @server.resource(
        "identifiers://{registry}",
        name="read_identifiers_registry",
        title="Identifiers.org registry namespace",
        description=(
            "Return one complete native identifiers.org namespace record, including "
            "all provider resources and institutions."
        ),
        mime_type="application/json",
    )
    async def read_registry(registry: str) -> str:
        """Return one official nested namespace object.

        Raises:
            ResourceError: If the registry prefix is unknown.
        """
        record = namespaces.get(registry.casefold())
        if record is None:
            raise ResourceError(f"Unknown identifiers.org registry {registry!r}")
        return _json_text(record)

    @server.tool(
        name="parse_identifiers",
        title="Parse identifiers.org IDs",
        description=(
            "Extract identifiers.org Compact Identifiers, identifiers.org URLs, "
            "BioV resource URIs, and allowlisted unambiguous bare IDs from a prompt; "
            "validate them and return MCP resource links."
        ),
        annotations=ToolAnnotations(
            read_only_hint=True,
            destructive_hint=False,
            idempotent_hint=True,
            open_world_hint=True,
        ),
    )
    async def parse_identifiers(prompt: str) -> list[TextContent | ResourceLink]:
        """Return resource links for resolver-valid identifiers in a prompt.

        Args:
            prompt: Complete prompt to scan for explicit identifiers.

        Returns:
            A summary followed by resource links, or a no-match message.

        Raises:
            ToolError: If identifiers.org cannot validate candidates.
        """
        links: list[ResourceLink] = []
        linked_uris: set[str] = set()
        for compact_id in extract_identifier_candidates(prompt):
            try:
                resolution = await resolve(compact_id)
                namespace_prefix, local_id = _parsed_identifier(resolution)
            except IdentifierNotFoundError:
                continue
            except IdentifierServiceError as error:
                raise ToolError(str(error)) from error
            namespace = namespaces.get(namespace_prefix.casefold())
            if namespace is None:
                raise ToolError(
                    f"Resolver namespace {namespace_prefix!r} is absent from the "
                    "packaged identifiers.org registry asset"
                )
            try:
                _validate_accession(namespace, local_id)
            except ResourceError as error:
                raise ToolError(str(error)) from error
            resource_uri = build_namespace_resource_uri(namespace, local_id)
            if resource_uri in linked_uris:
                continue
            linked_uris.add(resource_uri)
            data_backed = namespace["prefix"] in DATA_RESOURCE_NAMESPACE_PREFIXES
            links.append(
                ResourceLink(
                    name=compact_id,
                    title=compact_id,
                    uri=resource_uri,
                    description=(
                        f"Canonical {namespace['name']} data for {compact_id}"
                        if data_backed
                        else f"Identifiers.org resolution for {compact_id}"
                    ),
                    mime_type="application/json",
                    _meta={
                        "identifiersOrgUrl": build_identifiers_org_url(compact_id),
                        "namespacePrefix": namespace["prefix"],
                        "identifierPattern": namespace["pattern"],
                        "resourceKind": "data" if data_backed else "resolution",
                    },
                )
            )

        if not links:
            return [TextContent(text="No resolver-valid identifiers.org IDs found.")]
        summary = TextContent(
            text=f"Found {len(links)} resolver-valid identifiers.org ID(s)."
        )
        return [summary, *links]

    @server.tool(
        name="resolve_identifiers",
        title="Read an identifier resource",
        description=(
            "Read a refseq.gcf://, uniprot://, or identifiers:// resource URI and "
            "return its content as an embedded resource for tool-only MCP clients."
        ),
        annotations=ToolAnnotations(
            read_only_hint=True,
            destructive_hint=False,
            idempotent_hint=True,
            open_world_hint=True,
        ),
    )
    async def resolve_identifiers(uri: str) -> EmbeddedResource:
        """Read one BioV identifier resource through the tool interface.

        Args:
            uri: Exact URI returned by ``parse_identifiers`` or listed resources.

        Returns:
            The same content as ``resources/read``, embedded in the tool result.

        Raises:
            ToolError: If the URI matches no resource or returns no content.
        """
        contents = cast(
            "list[ReadResourceContents]",
            list(await server.read_resource(uri)),
        )
        if not contents:
            raise ToolError(f"Resource {uri!r} returned no content")
        item = contents[0]
        return EmbeddedResource(
            resource=TextResourceContents(
                uri=uri,
                mime_type=item.mime_type,
                text=cast("str", item.content),
                _meta=item.meta,
            )
        )

    return server


mcp = create_mcp_server()


def main() -> None:
    """Run the BioV MCP server over standard input/output."""
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()


__all__ = [
    "BARE_IDENTIFIER_NAMESPACE_ALLOWLIST",
    "create_mcp_server",
    "main",
    "mcp",
]
