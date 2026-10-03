"""Managed analysis MCP contracts, independent of a scientific environment."""

import json
import threading

import anyio
import pytest
from mcp.server.lowlevel.helper_types import ReadResourceContents
from mcp.server.mcpserver.exceptions import ResourceError, ToolError
from mcp.types import BlobResourceContents, CallToolResult, ReadResourceResult

import biov.mcp as mcp_module
from biov.analysis import AnalysisRequest
from biov.mcp import create_mcp_server


def _request() -> dict:
    """Return the smallest request using the declared Python entry point."""
    return {
        "name": "Protein properties",
        "environment": "python",
        "code": "print('analysis stdout')",
        "inputs": {},
        "outputs": {"properties.csv": "csv"},
    }


def _summary(status: str = "succeeded") -> dict:
    """Return a persisted run summary without provider or environment I/O."""
    return {
        "record": "file:///managed/run/record.json",
        "name": "Protein properties",
        "status": status,
        "stage": "program",
        "exit_code": 1 if status == "failed" else 0,
        "inputs": {},
        "outputs": {},
        "checks": [],
        "diagnostic": "failed program" if status == "failed" else None,
        "logs": {},
        "environment": {},
    }


def test_analysis_tools_expose_request_schema_and_execute_off_event_loop(monkeypatch):
    """Validate inputs through the SDK and offload synchronous analysis work."""
    calling_thread = threading.get_ident()
    summary = _summary()
    seen = []

    def run(request):
        assert isinstance(request, AnalysisRequest)
        assert threading.get_ident() != calling_thread
        seen.append(request)
        return summary

    monkeypatch.setattr(mcp_module, "execute_analysis", run)
    server = create_mcp_server()

    async def check():
        tools = {tool.name: tool for tool in await server.list_tools()}
        tool = tools["run_analysis"]
        schema = tool.input_schema
        request_schema = schema["$defs"]["AnalysisRequest"]
        assert request_schema["additionalProperties"] is False
        assert {"name", "environment", "code", "inputs", "outputs"} <= set(
            request_schema["properties"]
        )
        assert tool.annotations is not None
        assert tool.annotations.read_only_hint is False
        assert tool.annotations.idempotent_hint is False
        query_tool = tools["inspect_analysis"]
        assert query_tool.annotations is not None
        assert query_tool.annotations.read_only_hint is True

        result = await server.call_tool("run_analysis", {"request": _request()})
        assert isinstance(result, CallToolResult)
        assert result.structured_content == summary
        assert not result.is_error
        assert result.content == []
        with pytest.raises(ToolError):
            await server.call_tool(
                "run_analysis", {"request": {**_request(), "unexpected": True}}
            )

    anyio.run(check)
    assert len(seen) == 1


def test_real_stdio_rejects_remote_analysis_before_creating_results(tmp_path):
    """Report unsupported remote execution instead of silently running locally."""
    import sys

    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    config = tmp_path / "config.toml"
    config.write_text("")
    results = tmp_path / "results"
    parameters = StdioServerParameters(
        command=sys.executable,
        args=["-c", "from biov.mcp import main; main()"],
        env={
            "BIOV_CONFIG": str(config),
            "BIOV_EXECUTION_HOST": "unconfigured-test-host.invalid",
            "BIOV_ANALYSIS_ROOT": str(results),
        },
    )

    async def check():
        with anyio.fail_after(30):
            async with stdio_client(parameters) as (reader, writer):
                async with ClientSession(reader, writer) as session:
                    await session.initialize()
                    result = await session.call_tool(
                        "run_analysis", {"request": _request()}
                    )
                    assert isinstance(result, CallToolResult)
                    assert result.is_error
                    assert result.structured_content is not None
                    assert result.structured_content["stage"] == "launch"
                    assert result.structured_content["diagnostic"] == (
                        "Managed analysis currently supports local execution only"
                    )

    anyio.run(check)
    assert not results.exists()


def test_failed_run_and_successful_query_have_different_error_semantics(monkeypatch):
    """Report program failure without making a status query a tool error."""
    summary = _summary("failed")
    monkeypatch.setattr(mcp_module, "execute_analysis", lambda request: summary)
    monkeypatch.setattr(mcp_module, "inspect_saved_analysis", lambda record: summary)
    server = create_mcp_server()

    async def check():
        result = await server.call_tool("run_analysis", {"request": _request()})
        assert isinstance(result, CallToolResult)
        assert result.is_error
        query = await server.call_tool(
            "inspect_analysis", {"record": summary["record"]}
        )
        assert isinstance(query, CallToolResult)
        assert not query.is_error
        assert query.structured_content == result.structured_content == summary

    anyio.run(check)


@pytest.mark.parametrize(
    ("tool", "core", "arguments", "stage"),
    [
        ("run_analysis", "execute_analysis", {"request": _request()}, "launch"),
        (
            "inspect_analysis",
            "inspect_saved_analysis",
            {"record": "file:///missing/record.json"},
            "data_checks",
        ),
    ],
)
def test_analysis_boundary_errors_remain_bounded_and_structured(
    monkeypatch, tool, core, arguments, stage
):
    """Preserve expected access errors without exposing unbounded diagnostics."""

    def fail(value):
        raise FileNotFoundError("unavailable " + "x" * 40_000)

    monkeypatch.setattr(mcp_module, core, fail)
    server = create_mcp_server()

    async def call():
        return await server.call_tool(tool, arguments)

    result = anyio.run(call)
    assert isinstance(result, CallToolResult)
    assert result.is_error
    assert result.structured_content is not None
    assert result.structured_content["status"] == "failed"
    assert result.structured_content["stage"] == stage
    assert result.structured_content["diagnostic"].startswith("unavailable ")
    assert result.structured_content["diagnostic"].endswith("(truncated)")
    assert len(result.model_dump_json().encode()) < 32_768


def test_analysis_file_resource_preserves_bytes_and_core_access_checks(
    monkeypatch, tmp_path
):
    """Delegate managed-file authorization and limits without decoding bytes."""
    target = tmp_path / "protein data.bin"
    calls = []
    payload = b"\x00\xff\r\ncomplete output"

    def read(uri):
        calls.append(uri)
        if uri != target.as_uri():
            raise ValueError("Not a registered analysis file")
        return payload, "application/octet-stream"

    monkeypatch.setattr(mcp_module, "read_analysis_file", read)
    server = create_mcp_server()

    async def check():
        contents = list(await server.read_resource(target.as_uri()))
        content = contents[0]
        assert isinstance(content, ReadResourceContents)
        assert content.content == payload
        assert content.mime_type == "application/octet-stream"
        with pytest.raises(ResourceError, match="Not a registered analysis file"):
            await server.read_resource((tmp_path / "unregistered").as_uri())
        with pytest.raises(ResourceError):
            await server.read_resource(f"{tmp_path.as_uri()}/%2E%2E/private")

    anyio.run(check)
    assert calls == [
        target.as_uri(),
        (tmp_path / "unregistered").as_uri(),
        f"{tmp_path.as_uri()}/../private",
    ]


def test_analysis_stdio_serializes_structured_failure_and_file_bytes(tmp_path):
    """Use an actual SDK client and child server without claiming Codex support."""
    import base64
    import sys

    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client
    from mcp.shared.exceptions import MCPError

    script = tmp_path / "analysis_server.py"
    summary = _summary("failed")
    script.write_text(
        "import biov.mcp as m\n"
        f"summary = {summary!r}\n"
        "m.execute_analysis = lambda request: summary\n"
        "m.inspect_saved_analysis = lambda record: summary\n"
        "def read(uri):\n"
        "    if uri != 'file:///managed/run/result.bin':\n"
        "        raise ValueError('Not a registered analysis file')\n"
        "    return b'complete\\x00data', 'application/octet-stream'\n"
        "m.read_analysis_file = read\n"
        "m.create_mcp_server().run()\n"
    )

    async def check():
        with anyio.fail_after(30):
            async with stdio_client(
                StdioServerParameters(command=sys.executable, args=[str(script)])
            ) as (reader, writer):
                async with ClientSession(reader, writer) as session:
                    await session.initialize()
                    result = await session.call_tool(
                        "run_analysis", {"request": _request()}
                    )
                    assert isinstance(result, CallToolResult)
                    assert result.is_error
                    assert result.structured_content == summary
                    query = await session.call_tool(
                        "inspect_analysis", {"record": summary["record"]}
                    )
                    assert isinstance(query, CallToolResult)
                    assert not query.is_error
                    assert query.structured_content == summary
                    invalid = await session.call_tool(
                        "run_analysis", {"request": {**_request(), "unknown": 1}}
                    )
                    assert isinstance(invalid, CallToolResult)
                    assert invalid.is_error
                    resource = await session.read_resource(
                        "file:///managed/run/result.bin"
                    )
                    assert isinstance(resource, ReadResourceResult)
                    content = resource.contents[0]
                    assert isinstance(content, BlobResourceContents)
                    assert base64.b64decode(content.blob) == b"complete\0data"
                    with pytest.raises(MCPError, match="Not a registered"):
                        await session.read_resource("file:///private/secret")
                    assert len(json.dumps(result.model_dump(mode="json"))) < 32_768

    anyio.run(check)
