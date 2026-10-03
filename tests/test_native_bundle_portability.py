"""Opt-in installed-native export acceptance with an independent offline reader.

Configure both BIOV_TEST_BINARY and BIOV_STANDALONE_PYTHON to run this gate.
The latter must name a separate virtualenv containing only pyarrow==25.0.1,
with no BioV installation. CI provisions it separately from the pytest driver.
The consumer receives only the four moved bundle files, never an MCP response.
"""

import asyncio
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

_TIMEOUT = 30
_ROOT = Path(__file__).resolve().parents[1]
_CSV = """sample,text,count,ratio,active,drop_me
00123,000001,9007199254740993,1.25,true,excluded-column
00002,,,,false,excluded-column
00002,,,,false,excluded-column
00003,00003,20,2.5,,excluded-column
00004,9e3,-9,-0.0,true,excluded-column
00005,00100,-9007199254740993,-3.75,false,excluded-column
00006,00000,0,0.125,true,excluded-column
skip,not-exported,777,777,true,excluded-column
"""

# This program runs under the independent interpreter outside the checkout.
# It imports no BioV code and cannot call the producing server, SQLite, or a
# network service. Expected rows are declared here, not copied from its preview.
_STANDALONE_CHECK = r"""
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import sys

assert os.environ["PATH"] == ""
assert os.environ["PYTHONPATH"] == ""
assert sys.prefix != sys.base_prefix, "Reader must be a separate virtualenv"
assert importlib.util.find_spec("biov") is None, "BioV must not be installed"
assert {
    distribution.metadata["Name"].lower()
    for distribution in importlib.metadata.distributions()
} == {"pyarrow"}, "Provision a clean pyarrow-only reader virtualenv"

def offline_only(event, args):
    if event.startswith(("socket.", "subprocess.", "sqlite3.")) or event in {
        "os.system", "os.exec", "os.posix_spawn", "os.fork", "os.forkpty"
    }:
        raise AssertionError(f"Standalone reading must be offline: {event}")
    if event == "import" and args[0].split(".")[0] in {"biov", "sqlite3"}:
        raise AssertionError(f"Forbidden reader dependency: {args[0]}")

sys.addaudithook(offline_only)
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.ipc as ipc

assert pa.__version__ == "25.0.1", pa.__version__
root = Path.cwd()
manifests = list(root.glob("artifact_*.manifest.json"))
readmes = list(root.glob("artifact_*.README.md"))
assert len(manifests) == len(readmes) == 1, "Discover one conventional bundle"
manifest_path, readme_path = manifests[0], readmes[0]
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
assert manifest["manifest_version"] == 1
assert manifest["format"] == "arrow_ipc_file"
assert set(manifest["files"]) == {"arrow", "record", "readme"}

paths = {}
for role, filename in manifest["files"].items():
    relative = Path(filename)
    assert not relative.is_absolute() and relative.parts == (filename,)
    assert "\\" not in filename and filename not in {".", ".."}
    path = root / relative
    assert path.is_file() and not path.is_symlink(), (role, filename)
    paths[role] = path
assert paths["readme"] == readme_path
assert set(root.iterdir()) == {manifest_path, *paths.values()}, "Only four files"
stem = manifest_path.name.removesuffix(".manifest.json")
assert stem == manifest["artifact_id"]
assert paths["arrow"].name == stem + ".arrow"
assert paths["record"].name == stem + ".json"
assert readme_path.name == stem + ".README.md"

arrow_bytes = paths["arrow"].read_bytes()
record_bytes = paths["record"].read_bytes()
content = manifest["content"]
assert len(arrow_bytes) == content["bytes"]
assert hashlib.sha256(arrow_bytes).hexdigest() == content["sha256"]
assert len(record_bytes) == manifest["record"]["bytes"]
assert hashlib.sha256(record_bytes).hexdigest() == manifest["record"]["sha256"]
record = json.loads(record_bytes)
assert record["record_version"] == manifest["record"]["record_version"] == 2
assert set(record) == {
    "record_version", "artifact_id", "format", "file", "bytes", "sha256",
    "row_count", "schema", "scientific_metadata", "metadata_status",
    "provenance", "reopen_verification", "software"
}, "Companions must not extend the strict version-2 reopen record"
assert record["artifact_id"] == manifest["artifact_id"]
assert record["format"] == "arrow_ipc"
assert record["file"] == paths["arrow"].name
assert record["bytes"] == content["bytes"]
assert record["sha256"] == content["sha256"]
assert record["row_count"] == content["row_count"] == 7

# Actually run the shipped quick-start, rather than validating only a second
# implementation. Any forbidden dependency is rejected by the audit hook.
readme = readme_path.read_text(encoding="utf-8")
examples = re.findall(r"^```python\s*\n(.*?)^```\s*$", readme, re.M | re.S)
assert len(examples) == 1, "README must provide one executable Python quick-start"
exec(compile(examples[0], f"{readme_path} [Python quick-start]", "exec"),
     {"__name__": "__main__"})

# Read Arrow again independently of names/objects created by the README.
with pa.BufferReader(arrow_bytes) as reader:
    table = ipc.open_file(reader).read_all()
names = ["sample", "text", "count", "ratio", "active"]
logical_types = ["string", "string", "int64", "float64", "boolean"]
polars_types = ["str", "str", "i64", "f64", "bool"]
null_counts = [0, 2, 2, 2, 1]
assert table.column_names == names
assert table.num_rows == 7 and table.num_columns == 5
assert record["schema"] == [
    {"name": name, "dtype": dtype}
    for name, dtype in zip(names, polars_types, strict=True)
]
assert len(manifest["columns"]) == len(names)
for index, (field, column) in enumerate(zip(
    table.schema, manifest["columns"], strict=True
)):
    assert field.name == column["name"] == names[index]
    assert column["logical_type"] == logical_types[index]
    assert column["polars_dtype"] == polars_types[index]
    assert field.nullable and column["nullable"] is True
    assert table.column(index).null_count == column["null_count"] == null_counts[index]
    assert column["description"] is None
    assert column["units"] is None and column["coordinates"] is None
    if index < 2:
        assert field.type == pa.large_string(), "New exports must support Arrow kernels directly"
    else:
        assert field.type == (pa.int64(), pa.float64(), pa.bool_())[index - 2]

expected = [
    {"sample": "00002", "text": None, "count": None, "ratio": None, "active": False},
    {"sample": "00002", "text": None, "count": None, "ratio": None, "active": False},
    {"sample": "00003", "text": "00003", "count": 20, "ratio": 2.5, "active": None},
    {"sample": "00004", "text": "9e3", "count": -9, "ratio": -0.0, "active": True},
    {"sample": "00005", "text": "00100", "count": -9007199254740993, "ratio": -3.75, "active": False},
    {"sample": "00006", "text": "00000", "count": 0, "ratio": 0.125, "active": True},
    {"sample": "00123", "text": "000001", "count": 9007199254740993, "ratio": 1.25, "active": True},
]
assert table.to_pylist() == expected, table.to_pylist()
assert math.copysign(1, table["ratio"][3].as_py()) == -1, "Keep signed zero"
# New exports work with standard kernels immediately, without a compatibility
# cast, BioV helper, original CSV, or running server.
selected = table.filter(pc.greater(table["count"], pa.scalar(10, pa.int64())))
assert selected.to_pylist() == [expected[2], expected[6]]
assert pc.sum(table["count"]).as_py() == 11
assert pc.count(table["count"]).as_py() == 5
assert math.isclose(pc.mean(table["ratio"]).as_py(), 0.025, rel_tol=1e-14)
assert table.filter(pc.is_null(table["count"])).to_pylist() == expected[:2]
assert pc.sum(table["active"]).as_py() == 3

# Declared biological context is retained, without inferring unknown versions.
assert manifest["scientific_metadata"] == record["scientific_metadata"] == {
    "identifier": "refseq.gcf:GCF_000001405.40",
    "species": None,
    "reference": "synthetic fixture reference",
    "coordinates": None,
    "units": "synthetic count",
}
identifier = manifest["biological_identifier"]
assert identifier["namespace"] == "refseq.gcf"
assert identifier["accession"] == "GCF_000001405.40"
assert identifier["base_accession"] == "GCF_000001405"
assert str(identifier["accession_version"]) == "40"
assert identifier["validation"] == "syntax_only_not_provider_verified"
for key in ("entry_version", "sequence_version", "provider_release"):
    assert identifier[key] is None
assert manifest["versions"] == {"reference_version": None, "provider_release": None}
source = manifest["source"]
assert source["historical_path"] == "typed.csv"
assert source["path_interpretation"] == "relative_to_original_data_root_not_bundle"
assert source["required_for_reading"] is False
assert source["sha256"] == record["provenance"]["source_sha256"]
assert source["bytes"] == record["provenance"]["source_bytes"]
assert not (root / source["historical_path"]).exists()
lineage = manifest["lineage"]
assert lineage["historical_references_only"] is True
assert lineage["operations"] == record["provenance"]["operations"]
assert len(lineage["operations"]) == 1
assert lineage["declared_schema"] == {"count": "int64", "ratio": "float64", "active": "boolean"}
assert lineage["reopen_verification"] is None
assert manifest["software"] == record["software"]
assert not any(name == "biov" or name.startswith("biov.") for name in sys.modules)
assert "sqlite3" not in sys.modules
print("standalone moved-bundle filter and summary passed")
"""


async def _export_bundle(binary: Path, data: Path, output: Path) -> dict:
    """Run the installed MCP server and stop it before returning export metadata.

    Returns:
        Structured native export result after successful server shutdown.
    """
    environment = {"PATH": "", "PYTHONPATH": "", "POLARS_MAX_THREADS": "2"}
    process = await asyncio.wait_for(
        asyncio.create_subprocess_exec(
            str(binary),
            "mcp",
            "--data-root",
            str(data),
            "--output-root",
            str(output),
            cwd=data.parent,
            env=environment,
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            limit=128 * 1024,
        ),
        timeout=_TIMEOUT,
    )
    assert process.stdin is not None and process.stdout is not None
    sequence = 0

    async def request(method: str, params: dict) -> dict:
        """Return a matching JSON-RPC result within a bounded request deadline."""
        nonlocal sequence
        sequence += 1
        message = {"jsonrpc": "2.0", "id": sequence, "method": method, "params": params}
        assert process.stdin is not None and process.stdout is not None
        process.stdin.write(json.dumps(message).encode() + b"\n")
        await process.stdin.drain()
        while True:
            line = await process.stdout.readline()
            assert line, f"MCP exited before replying to {method}"
            response = json.loads(line)
            assert response["jsonrpc"] == "2.0", response
            if "id" not in response:
                assert "method" in response, response
                continue
            assert response["id"] == sequence and "error" not in response, response
            return response["result"]

    async def call(name: str, arguments: dict) -> dict:
        """Validate an actual native tool response, without importing BioV.

        Returns:
            The tool's structured result.
        """
        result = await asyncio.wait_for(
            request("tools/call", {"name": name, "arguments": arguments}),
            timeout=_TIMEOUT,
        )
        assert result.get("isError") is False, (name, result)
        assert isinstance(result.get("structuredContent"), dict), result
        return result["structuredContent"]

    try:
        initialized = await asyncio.wait_for(
            request(
                "initialize",
                {
                    "protocolVersion": "2025-06-18",
                    "capabilities": {},
                    "clientInfo": {
                        "name": "independent-bundle-acceptance",
                        "version": "1",
                    },
                },
            ),
            timeout=_TIMEOUT,
        )
        assert initialized["serverInfo"]["name"] == "biov-rs"
        process.stdin.write(b'{"jsonrpc":"2.0","method":"notifications/initialized"}\n')
        await asyncio.wait_for(process.stdin.drain(), timeout=_TIMEOUT)
        opened = await call(
            "dataset_open",
            {
                "path": "typed.csv",
                "schema": {"count": "int64", "ratio": "float64", "active": "boolean"},
                "metadata": {
                    "identifier": "refseq.gcf:GCF_000001405.40",
                    "reference": "synthetic fixture reference",
                    "units": "synthetic count",
                },
                "preview_rows": 1,
            },
        )
        assert opened["row_count"] == 8
        derived = await call(
            "dataset_query",
            {
                "dataset_id": opened["dataset_id"],
                "filter": {"column": "sample", "op": "lt", "value": "01000"},
                "sort": {"column": "sample", "descending": False},
                "select": ["sample", "text", "count", "ratio", "active"],
                "preview_rows": 1,
            },
        )
        assert derived["dataset_id"] != opened["dataset_id"]
        assert derived["row_count"] == 7
        assert derived["preview"]["returned_rows"] == 1
        exported = await call("dataset_export", {"dataset_id": derived["dataset_id"]})
        process.stdin.close()
        stdout, stderr = await asyncio.wait_for(process.communicate(), timeout=_TIMEOUT)
        assert process.returncode == 0, stderr.decode(errors="replace")
        assert not stdout and not stderr, (stdout, stderr)
        return exported
    finally:
        if process.returncode is None:
            process.kill()
            await asyncio.wait_for(process.communicate(), timeout=_TIMEOUT)


def _configured_executable(name: str) -> Path:
    """Reject explicit broken acceptance configuration instead of silently skipping.

    Returns:
        The configured executable's absolute path without resolving venv symlinks.
    """
    configured = os.environ.get(name)
    assert configured, (
        f"Set both BIOV_TEST_BINARY and BIOV_STANDALONE_PYTHON; missing {name}"
    )
    # Keep the virtualenv executable path; resolving its symlink loses the venv.
    path = Path(configured).absolute()
    assert path.is_file() and os.access(path, os.X_OK), (
        f"{name} is not executable: {path}"
    )
    return path


def test_native_bundle_survives_offline_relocation(tmp_path: Path) -> None:
    """Analyze a moved native bundle using only its files and stock pinned PyArrow."""
    if not any(
        name in os.environ for name in ("BIOV_TEST_BINARY", "BIOV_STANDALONE_PYTHON")
    ):
        pytest.skip(
            "Set BIOV_TEST_BINARY and BIOV_STANDALONE_PYTHON for installed acceptance"
        )
    binary = _configured_executable("BIOV_TEST_BINARY")
    python = _configured_executable("BIOV_STANDALONE_PYTHON")
    data = tmp_path / "original-data"
    output = tmp_path / "original-output"
    data.mkdir()
    output.mkdir()
    (data / "typed.csv").write_text(_CSV, encoding="utf-8")
    exported = asyncio.run(_export_bundle(binary, data, output))
    source_identity = exported["record"]["provenance"]
    assert source_identity["source_bytes"] == len(_CSV.encode())
    assert source_identity["source_sha256"] == hashlib.sha256(_CSV.encode()).hexdigest()
    keys = ("execution_host_path", "record_path", "manifest_path", "readme_path")
    assert all(key in exported for key in keys), (
        "Native export must return four bundle paths"
    )
    artifacts = [Path(exported[key]) for key in keys]
    assert all(path.is_file() and path.parent == output for path in artifacts)
    assert set(output.iterdir()) == set(artifacts), (
        "Export must leave exactly four files"
    )

    # A different temporary root prevents accidental access through the producer
    # working directory. No source, helper, side database, or response is copied.
    with tempfile.TemporaryDirectory(prefix="independent-arrow-consumer-") as directory:
        consumer = Path(directory).resolve()
        assert _ROOT != consumer and _ROOT not in consumer.parents
        assert tmp_path != consumer and tmp_path not in consumer.parents
        for artifact in artifacts:
            shutil.move(artifact, consumer / artifact.name)
        shutil.rmtree(data)
        shutil.rmtree(output)
        assert not data.exists() and not output.exists()
        try:
            completed = subprocess.run(
                [str(python), "-I", "-c", _STANDALONE_CHECK],
                cwd=consumer,
                env={"PATH": "", "PYTHONPATH": "", "PYTHONNOUSERSITE": "1"},
                capture_output=True,
                text=True,
                timeout=60,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            pytest.fail(f"Standalone PyArrow consumer exceeded 60 seconds: {error}")
        assert completed.returncode == 0, (
            "Independent moved-bundle analysis failed. Provision the reader with only "
            "pyarrow==25.0.1 and rebuild the installed native binary.\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
        assert "standalone moved-bundle filter and summary passed" in completed.stdout
