"""BioV Model Context Protocol server."""

import json
import re
import subprocess  # noqa: S404 - only exception types
from collections.abc import Awaitable, Callable
from functools import partial, wraps
from pathlib import Path
from typing import cast

from anyio.to_thread import run_sync
from mcp.server import MCPServer
from mcp.server.lowlevel.helper_types import ReadResourceContents
from mcp.server.mcpserver.exceptions import ResourceError, ToolError
from mcp.server.mcpserver.resources.templates import ResourceSecurity
from mcp.types import (
    CallToolResult,
    EmbeddedResource,
    ResourceLink,
    TextContent,
    TextResourceContents,
    ToolAnnotations,
)

from .analysis import (
    AnalysisRequest,
    read_analysis_file,
)
from .analysis import (
    inspect_analysis as inspect_saved_analysis,
)
from .analysis import (
    run_analysis as execute_analysis,
)
from .artifacts import (
    Artifact,
    artifact_capabilities,
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


def _resource_errors[**P, T](
    handler: Callable[P, Awaitable[T]],
) -> Callable[P, Awaitable[T]]:
    """Expose expected provider failures at the MCP boundary only.

    Returns:
        Handler using the SDK error type for expected provider failures.
    """

    @wraps(handler)
    async def call(*args: P.args, **kwargs: P.kwargs) -> T:
        try:
            return await handler(*args, **kwargs)
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
            raise ResourceError(str(error)) from error

    return call


def _json_text(value: object) -> str:
    """Serialize one native upstream JSON value for MCP text transport.

    Returns:
        Indented UTF-8 JSON text without changing its data model.
    """
    return json.dumps(value, ensure_ascii=False, indent=2)


def _analysis_error(error: OSError | ValueError, stage: str) -> CallToolResult:
    """Describe an expected analysis boundary error without hiding its status.

    Returns:
        A bounded diagnostic with explicit MCP tool-error semantics.
    """
    diagnostic = str(error)
    if len(diagnostic) > 4096:
        diagnostic = diagnostic[:4096] + "… (truncated)"
    return CallToolResult(
        content=[],
        structured_content={
            "status": "failed",
            "stage": stage,
            "diagnostic": diagnostic,
        },
        is_error=True,
    )


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

    @_resource_errors
    async def read_data(accession: str) -> str:
        """Return the provider-native canonical JSON representation."""
        _validate_accession(namespace, accession)
        return await run_sync(reader, accession)

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
        Configured server with data descriptions, identifier resources, and tools.
    """
    resolve = resolver or resolve_identifier
    resource_path = artifact_path or path
    read_refseq = refseq_reader or genome_summary
    namespaces = namespaces_by_prefix()
    server = MCPServer(
        name="biov",
        title="BioV",
        description="Biological data and managed analysis for AI agents",
        instructions=(
            "Call parse_identifiers with prompts containing biological IDs. "
            "Read refseq.gcf://, uniprot://, pubmed://, clinvar://, dbsnp://, "
            "and geo:// resources for upstream data. Read identifiers://<registry> "
            "for native registry metadata and "
            "identifiers://<registry>:<id> for native resolver JSON. Tool-only "
            "clients can pass identifier "
            "resource URIs to resolve_identifiers. "
            "For RefSeq and UniProt files, generated Python should pass "
            "biov.path(uri) to libraries accepting "
            "os.PathLike, or use biov.open(uri, mode='rb' or 'rt') when a readable "
            "file object is required. Call both only inside the execution environment. "
            "Other identifier resources describe downloadable analysis files and "
            "their formats. Use biov.path, biov.open, or fsspec.open on their URI "
            "inside the execution environment to obtain the file. "
            "Call run_analysis with a complete ordinary Python script, named "
            "inputs, declared outputs and a configured locked Pixi environment "
            "to execute an analysis directly through MCP. Reuse complete output "
            "file references as subsequent inputs, never preview rows. Call "
            "inspect_analysis with the returned record to inspect saved facts; "
            "unknown status must not trigger automatic resubmission. Managed "
            "file resources return at most 1 MiB; use an output's download_url "
            "for larger files and verify its byte size and SHA-256. File URIs "
            "refer to this server's storage, not the client's filesystem. "
            "The BioV plugin includes biological-data, scientific-software, and "
            "task-specific skills; load their instructions and references on demand "
            "through the host's skill mechanism. Use scientific-software for the "
            "software environment and execution instructions. Run scripts or native "
            "commands with biov exec COMMAND ARGS... through the host's terminal "
            "tool. Execution environments are configured separately; a listed "
            "software package is not necessarily installed."
        ),
    )

    def read_uniprot_json(accession: str) -> str:
        """Read the original cached UniProt entry JSON.

        Returns:
            Original provider JSON text.
        """
        artifact = resource_path(f"uniprot://{accession}", artifact="entry_json")
        return artifact.path.read_text(encoding="utf-8")

    readers = [
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
    ]

    def describe_file(accession: str, prefix: str) -> str:
        """Return available representations for one locally validated file ID."""
        capability = artifact_capabilities()["namespaces"][prefix]
        return _json_text(
            {
                "uri": build_namespace_resource_uri(namespaces[prefix], accession),
                "default_kind": capability["default_kind"],
                "representations": {
                    kind: {
                        key: value for key, value in properties.items() if key != "url"
                    }
                    for kind, properties in capability["kinds"].items()
                },
            }
        )

    readers.extend(
        (
            prefix,
            "Describe the analysis file; open its URI through BioV or fsspec.",
            partial(describe_file, prefix=prefix),
        )
        for prefix in sorted(
            DATA_RESOURCE_NAMESPACE_PREFIXES - {"refseq.gcf", "uniprot"}
        )
    )
    for prefix, description, reader in readers:
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
    @_resource_errors
    async def read_identifier(registry: str, id: str) -> str:
        """Return the official resolver JSON for one registry-local ID.

        Raises:
            ResourceError: If registry or accession validation fails.
        """
        namespace = namespaces.get(registry.casefold())
        if namespace is None:
            raise ResourceError(f"Unknown identifiers.org registry {registry!r}")
        _validate_accession(namespace, id)
        resolution = await resolve(_compact_identifier(namespace, id))
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
            "validate them and return a JSON summary with MCP resource links."
        ),
        annotations=ToolAnnotations(
            read_only_hint=True,
            destructive_hint=False,
            idempotent_hint=True,
            open_world_hint=True,
        ),
    )
    @_resource_errors
    async def parse_identifiers(prompt: str) -> list[TextContent | ResourceLink]:
        """Return resource links for resolver-valid identifiers in a prompt.

        Args:
            prompt: Complete prompt to scan for explicit identifiers.

        Returns:
            A JSON summary followed by resource links.

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
            namespace = namespaces.get(namespace_prefix.casefold())
            if namespace is None:
                raise ToolError(
                    f"Resolver namespace {namespace_prefix!r} is absent from the "
                    "packaged identifiers.org registry asset"
                )
            _validate_accession(namespace, local_id)
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

        summary = TextContent(
            text=_json_text(
                {
                    "count": len(links),
                    "resource_uris": [str(link.uri) for link in links],
                }
            )
        )
        return [summary, *links]

    @server.tool(
        name="resolve_identifiers",
        title="Read an identifier resource",
        description=(
            "Read an identifier resource URI and "
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

    @server.tool(
        name="run_analysis",
        title="Run a managed analysis",
        description=(
            "Run a complete Python script in a declared locked Pixi environment "
            "on this server. Save complete outputs and execution records, and "
            "return bounded previews with reusable file references. The script "
            "receives inputs.json and parameters.json as arguments."
        ),
        annotations=ToolAnnotations(
            read_only_hint=False,
            destructive_hint=True,
            idempotent_hint=False,
            open_world_hint=True,
        ),
    )
    async def run_analysis(request: AnalysisRequest) -> CallToolResult:
        """Execute one analysis and preserve tool-error semantics for failure.

        Returns:
            Bounded structured result, with failed runs marked as tool errors.
        """
        try:
            summary = await run_sync(execute_analysis, request)
        except (OSError, ValueError) as error:
            return _analysis_error(error, "launch")
        return CallToolResult(
            content=[],
            structured_content=summary,
            is_error=summary["status"] == "failed",
        )

    @server.tool(
        name="inspect_analysis",
        title="Inspect a saved analysis",
        description=(
            "Read a saved analysis record and bounded previews without executing "
            "or resubmitting it. A failed analysis is a successful status query."
        ),
        annotations=ToolAnnotations(
            read_only_hint=True,
            destructive_hint=False,
            idempotent_hint=True,
            open_world_hint=False,
        ),
    )
    async def inspect_analysis(record: str) -> CallToolResult:
        """Read persisted facts without treating a failed run as a failed query.

        Returns:
            Bounded structured result or an explicit record-access error.
        """
        try:
            summary = await run_sync(inspect_saved_analysis, record)
        except (OSError, ValueError) as error:
            return _analysis_error(error, "data_checks")
        return CallToolResult(content=[], structured_content=summary)

    @server.resource(
        "file://{+path}",
        name="read_analysis_file",
        title="Saved analysis file",
        description=(
            "Read a registered analysis output, record or log of at most 1 MiB. "
            "Only files in the configured analysis root are accessible."
        ),
        mime_type="application/octet-stream",
        security=ResourceSecurity(reject_absolute_paths=False),
    )
    @_resource_errors
    async def read_saved_file(path: str) -> bytes:
        """Return the complete bytes of one validated managed file.

        Returns:
            Original file bytes after the core's identity and size checks.
        """
        content, _mime_type = await run_sync(read_analysis_file, Path(path).as_uri())
        return content

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
