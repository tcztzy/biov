"""Opt-in real Codex client verification, distinct from SDK protocol tests."""

import functools
import hashlib
import http.server
import json
import os
import shutil
import subprocess
import sys
import threading
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "docs/examples/sequence-analysis"


@pytest.mark.skipif(
    os.environ.get("BIOV_TEST_CODEX") != "1",
    reason="Set BIOV_TEST_CODEX=1 to use the authenticated Codex CLI and real Pixi",
)
def test_codex_calls_analysis_and_downloads_complete_results(tmp_path: Path) -> None:
    """Make the real client execute, query, and retrieve results over HTTP."""
    executable = shutil.which("codex")
    assert executable is not None, "Codex CLI is required for this explicit test"
    results = tmp_path / "server-results"
    results.mkdir()
    downloads = tmp_path / "client-downloads"
    downloads.mkdir()
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=results)
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base_url = f"http://127.0.0.1:{server.server_port}/"
    config = {
        "command": sys.executable,
        "args": ["-c", "from biov.mcp import main; main()"],
        "env": {
            "PATH": os.environ["PATH"],
            "BIOV_ENVIRONMENT_MANIFEST": str(EXAMPLE / "pyproject.toml"),
            "BIOV_ANALYSIS_ROOT": str(results),
            "BIOV_ANALYSIS_BASE_URL": base_url,
        },
        "startup_timeout_sec": 30,
        "tool_timeout_sec": 180,
    }
    # Read only server names: do not change persistent client settings or log auth.
    config_path = (
        Path(os.environ.get("CODEX_HOME", "~/.codex")).expanduser() / "config.toml"
    )
    configured = tomllib.loads(config_path.read_text()) if config_path.exists() else {}
    argv = [
        executable,
        "exec",
        "--ephemeral",
        "--json",
        "--approve-for-me",
        "-C",
        str(ROOT),
    ]
    for name in configured.get("mcp_servers", {}):
        argv.extend(["-c", f"mcp_servers.{name}.enabled=false"])
    for key, value in config.items():
        if isinstance(value, dict):
            for env_key, env_value in value.items():
                argv.extend(
                    [
                        "-c",
                        f"mcp_servers.biov_acceptance.env.{env_key}={json.dumps(env_value)}",
                    ]
                )
        else:
            argv.extend(
                ["-c", f"mcp_servers.biov_acceptance.{key}={json.dumps(value)}"]
            )
    first = {
        "name": "Translate five complete CDSs",
        "environment": "python",
        "code": (EXAMPLE / "extract_cds.py").read_text(),
        "inputs": {"genbank": str(ROOT / "tests/data/cor6_6.gb")},
        "outputs": {"proteins.fasta": "fasta", "cds.csv": "csv"},
        "parameters": {
            "accessions": ["X55053.1", "X62281.1", "M81224.1", "L31939.1", "AF297471.1"]
        },
        "requirements": [
            "Complete standard-code CDS and agreement with annotated translation"
        ],
    }
    second = {
        "name": "Compute properties for all five proteins",
        "environment": "python",
        "code": (EXAMPLE / "protein_properties.py").read_text(),
        "inputs": {"proteins": "REPLACE_WITH_FIRST_RESULT_URI"},
        "outputs": {"properties.csv": "csv"},
        "parameters": {},
        "requirements": ["Nonempty canonical amino acids with distinct identifiers"],
    }
    transfer = {
        "name": "Synthetic file transport check",
        "environment": "python",
        "code": "from pathlib import Path\nPath('transport.bin').write_bytes(b'BioV transport check\\n' * 60000)\n",
        "outputs": {"transport.bin": "file"},
    }
    prompt = (
        "You are testing this project's real MCP integration. Do not edit project files, "
        "commit, start other agents, or change settings. Use the biov_acceptance MCP "
        "run_analysis tool (not shell execution of analysis scripts) for each of these "
        "three exact requests in order. For the second request, use the complete "
        "proteins.fasta uri returned by the first in place of its input placeholder. "
        "Call inspect_analysis on the first returned record. Inspect the bounded "
        "previews. Then use a terminal Python urllib HTTP download to retrieve "
        "properties.csv and transport.bin from their returned download_url values "
        f"into {downloads}. Do not open files or paths in the server's results directory. "
        "All downloads must use HTTP even though this test runs on one machine. "
        "Verify each downloaded SHA-256 and byte size against the returned metadata; "
        "count the CSV's data rows (should be five). The synthetic transport fixture "
        "is deliberately larger than 1 MiB. Report observed statuses, row count and "
        "checksum comparisons concisely. Do not claim success if any step is blocked.\n"
        + json.dumps([first, second, transfer])
    )
    try:
        completed = subprocess.run(
            [*argv, "-"],
            input=prompt,
            text=True,
            capture_output=True,
            timeout=600,
            check=False,
        )
        (tmp_path / "codex.jsonl").write_text(completed.stdout)
        (tmp_path / "codex.stderr").write_text(completed.stderr)
        assert completed.returncode == 0, completed.stderr[-2000:]
        events = [
            json.loads(line)
            for line in completed.stdout.splitlines()
            if line.startswith("{")
        ]
        calls = [
            event.get("item", {})
            for event in events
            if event.get("type") == "item.completed"
        ]
        assert sum(item.get("tool") == "run_analysis" for item in calls) >= 3, (
            completed.stdout[-4000:]
        )
        assert any(item.get("tool") == "inspect_analysis" for item in calls), (
            completed.stdout[-4000:]
        )
        expected = b"BioV transport check\n" * 60000
        assert (downloads / "transport.bin").read_bytes() == expected
        properties = (downloads / "properties.csv").read_bytes()
        assert len(properties.decode().splitlines()) == 6
        records = [
            json.loads(path.read_text()) for path in results.glob("*/record.json")
        ]
        assert len(records) == 3
        assert all(record["status"] == "succeeded" for record in records)
        matching = [
            record for record in records if "properties.csv" in record["outputs"]
        ]
        assert (
            matching[0]["outputs"]["properties.csv"]["sha256"]
            == hashlib.sha256(properties).hexdigest()
        )
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
