"""Opt-in real Pixi, MCP stdio, and HTTP acceptance for the sequence example."""

import base64
import csv
import hashlib
import io
import json
import os
import sys
import threading
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlsplit
from urllib.request import urlopen

import anyio
import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.shared.exceptions import MCPError
from mcp.types import BlobResourceContents, CallToolResult, ReadResourceResult

pytestmark = pytest.mark.skipif(
    os.environ.get("BIOV_TEST_ANALYSIS") != "1",
    reason="Set BIOV_TEST_ANALYSIS=1 to run the installed real Pixi environment.",
)
ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "docs/examples/sequence-analysis"
ACCESSIONS = ["X55053.1", "X62281.1", "M81224.1", "L31939.1", "AF297471.1"]
PROTEINS = ["CAA38894.1", "CAA44171.1", "AAA32993.1", "AAA91051.1", "AAG13407.1"]


@asynccontextmanager
async def connect(parameters: StdioServerParameters) -> AsyncIterator[ClientSession]:
    """Start the actual BioV stdio server and establish an official SDK session.

    Yields:
        Initialized MCP client connected to a fresh server process.
    """
    async with (
        stdio_client(parameters) as (reader, writer),
        ClientSession(reader, writer, read_timeout_seconds=90) as session,
    ):
        await session.initialize()
        yield session


async def call(
    session: ClientSession, tool: str, arguments: dict, *, error: bool = False
) -> dict:
    """Check the real protocol's tool-error flag and bounded structured response.

    Returns:
        The tool's structured result.
    """
    result = await session.call_tool(tool, arguments)
    assert isinstance(result, CallToolResult)
    assert bool(result.is_error) is error, result
    assert len(result.model_dump_json().encode()) <= 32_768
    assert result.structured_content is not None
    return result.structured_content


async def read_bytes(session: ClientSession, uri: str) -> bytes:
    """Retrieve complete binary bytes through MCP, without using local file paths.

    Returns:
        Decoded resource bytes.
    """
    resource = await session.read_resource(uri)
    assert isinstance(resource, ReadResourceResult)
    assert len(resource.contents) == 1
    content = resource.contents[0]
    assert isinstance(content, BlobResourceContents)
    return base64.b64decode(content.blob, validate=True)


def download(url: str, destination: Path) -> bytes:
    """Download via the configured HTTP endpoint into independent client storage.

    Returns:
        Bytes read back from the downloaded client file.
    """
    with urlopen(url, timeout=10) as response:  # noqa: S310 - test-owned loopback URL
        destination.write_bytes(response.read())
    return destination.read_bytes()


def test_real_sequence_analysis_and_result_retrieval(tmp_path: Path) -> None:
    """Run real science, restart the server, transfer large results, and reject damage."""
    results_root = tmp_path / "results"
    results_root.mkdir()
    client_downloads = tmp_path / "client-downloads"
    client_downloads.mkdir()
    config = tmp_path / "config.toml"
    config.write_text("")
    server = ThreadingHTTPServer(
        ("127.0.0.1", 0),
        partial(SimpleHTTPRequestHandler, directory=str(results_root)),
    )
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base_url = f"http://127.0.0.1:{server.server_port}"
    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("BIOV_")
    }
    environment.update(
        BIOV_CONFIG=str(config),
        BIOV_ENVIRONMENT_MANIFEST=str(EXAMPLE / "pyproject.toml"),
        BIOV_ANALYSIS_ROOT=str(results_root),
        BIOV_ANALYSIS_BASE_URL=base_url,
    )
    parameters = StdioServerParameters(
        command=sys.executable,
        args=["-c", "from biov.mcp import main; main()"],
        env=environment,
        cwd=ROOT,
    )
    properties_code = (EXAMPLE / "protein_properties.py").read_text()
    extraction_code = (
        EXAMPLE / "extract_cds.py"
    ).read_text() + '\nprint("analysis stdout remains in the execution log")\n'
    first_request = {
        "name": "Translate five complete annotated CDSs",
        "environment": "python",
        "code": extraction_code,
        "inputs": {"genbank": str(ROOT / "tests/data/cor6_6.gb")},
        "outputs": {"proteins.fasta": "fasta", "cds.csv": "csv"},
        "parameters": {"accessions": ACCESSIONS},
        "requirements": [
            "Complete standard-code CDSs with matching annotated translation",
            "The researcher explicitly selects the five versioned accessions.",
        ],
    }

    async def exercise() -> None:
        async with connect(parameters) as session:
            names = {tool.name for tool in (await session.list_tools()).tools}
            assert {"run_analysis", "inspect_analysis"} <= names
            first = await call(session, "run_analysis", {"request": first_request})
            assert first["status"] == "succeeded"
            assert first["exit_code"] == 0
            assert first["checks"] and all(item["passed"] for item in first["checks"])
            protein_output = first["outputs"]["proteins.fasta"]
            assert protein_output["preview"]["truncated"] is True
            assert protein_output["preview"]["total_records"] is None
            assert [
                item["id"] for item in protein_output["preview"]["records"]
            ] == PROTEINS[:2]
            protein_bytes = await read_bytes(session, protein_output["uri"])
            assert hashlib.sha256(protein_bytes).hexdigest() == protein_output["sha256"]
            assert protein_bytes.count(b">") == 5
            first_record = await read_bytes(session, first["record"])
            recorded = json.loads(first_record)
            assert (
                recorded["code"]["sha256"]
                == hashlib.sha256(extraction_code.encode()).hexdigest()
            )
            assert recorded["parameters"] == {"accessions": ACCESSIONS}
            runtime = recorded["environment"]
            assert {
                "manifest_sha256",
                "lock_sha256",
                "pixi_version",
                "execution",
            } <= runtime.keys()
            assert runtime["execution"]["status"] == "completed"
            assert runtime["execution"]["exit_code"] == 0
            assert {
                name.lower(): value
                for name, value in runtime["execution"]["packages"].items()
            }["biopython"] == "1.88"
            assert b"analysis stdout remains in the execution log" in await read_bytes(
                session, first["logs"]["stdout"]
            )
            second_request = {
                "name": "Calculate all five proteins' properties",
                "environment": "python",
                "code": properties_code,
                "inputs": {"proteins": protein_output["uri"]},
                "outputs": {"properties.csv": "csv"},
                "parameters": {},
                "requirements": ["Average molecular weight and theoretical pI"],
            }
            second = await call(session, "run_analysis", {"request": second_request})
            assert second["status"] == "succeeded"
            assert second["record"] != first["record"]
            properties_output = second["outputs"]["properties.csv"]
            assert properties_output["preview"]["truncated"] is True
            assert properties_output["preview"]["total_rows"] is None
            assert len(properties_output["preview"]["rows"]) == 2
            property_bytes = await read_bytes(session, properties_output["uri"])
            assert (
                hashlib.sha256(property_bytes).hexdigest()
                == properties_output["sha256"]
            )
            rows = list(csv.DictReader(io.StringIO(property_bytes.decode())))
            assert [row["protein_id"] for row in rows] == PROTEINS
            assert [int(row["length_aa"]) for row in rows] == [66, 67, 65, 65, 65]
            assert [float(row["molecular_weight_da"]) for row in rows] == pytest.approx(
                [6551.1852, 7407.5077, 6552.2611, 6604.4217, 6536.2617]
            )
            assert [float(row["theoretical_pi"]) for row in rows] == pytest.approx(
                [9.1007600784, 10.9675710678, 9.1576856613, 9.1576856613, 9.1576856613]
            )
            failed = await call(
                session,
                "run_analysis",
                {
                    "request": {
                        **first_request,
                        "parameters": {"accessions": ["AJ237582.1"]},
                    }
                },
                error=True,
            )
            assert failed["status"] == "failed"
            assert failed["stage"] == "data_checks"
            assert any(not item["passed"] for item in failed["checks"])
            assert not failed["outputs"]
            # The next step also fails without touching the checked first result.
            failed_next = await call(
                session,
                "run_analysis",
                {"request": {**second_request, "parameters": {"unsupported": True}}},
                error=True,
            )
            assert failed_next["status"] == "failed"
            assert await read_bytes(session, first["record"]) == first_record
            assert await read_bytes(session, protein_output["uri"]) == protein_bytes

        # A new server resolves terminal facts without rerunning the programs.
        async with connect(parameters) as session:
            restored = await call(
                session, "inspect_analysis", {"record": first["record"]}
            )
            assert restored["status"] == "succeeded"
            queried_failure = await call(
                session, "inspect_analysis", {"record": failed["record"]}
            )
            assert queried_failure["status"] == "failed"
            assert await read_bytes(session, first["record"]) == first_record
            large = await call(
                session,
                "run_analysis",
                {
                    "request": {
                        "name": "Synthetic transport fixture, not biological data",
                        "environment": "python",
                        "code": "from pathlib import Path\nPath('payload.bin').write_bytes(bytes(range(256)) * 4097)\nprint('captured stdout')\n",
                        "inputs": {},
                        "outputs": {"payload.bin": "file"},
                    }
                },
            )
            output = large["outputs"]["payload.bin"]
            assert output["size"] > 1_048_576
            with pytest.raises(MCPError):
                await session.read_resource(output["uri"])
            assert output["download_url"].startswith(base_url + "/")
            destination = client_downloads / "payload.bin"
            downloaded = await anyio.to_thread.run_sync(
                download, output["download_url"], destination
            )
            assert len(downloaded) == output["size"]
            assert hashlib.sha256(downloaded).hexdigest() == output["sha256"]
            assert downloaded == bytes(range(256)) * 4097

            # Alter only this temporary test's outputs; nothing may repair them silently.
            protein_path = anyio.Path(unquote(urlsplit(protein_output["uri"]).path))
            await protein_path.write_bytes(b">changed\nMA\n")
            damaged = await call(
                session, "inspect_analysis", {"record": first["record"]}, error=True
            )
            assert damaged["status"] == "failed"
            reused = await call(
                session, "run_analysis", {"request": second_request}, error=True
            )
            assert reused["status"] == "failed"
            assert await protein_path.read_bytes() == b">changed\nMA\n"
            await protein_path.unlink()
            missing = await call(
                session, "inspect_analysis", {"record": first["record"]}, error=True
            )
            assert missing["status"] == "failed"
            missing_reuse = await call(
                session, "run_analysis", {"request": second_request}, error=True
            )
            assert missing_reuse["status"] == "failed"
            assert not await protein_path.exists()
            assert await read_bytes(session, first["record"]) == first_record
            assert await read_bytes(session, second["record"])

    try:
        anyio.run(exercise)
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
