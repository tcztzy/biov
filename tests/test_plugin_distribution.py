"""Check shared plugin contents and the declared MCP process entry point."""

import json
import os
import re
import sys
import tomllib
from pathlib import Path

import anyio
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.types import TextResourceContents

ROOT = Path(__file__).resolve().parents[1]


def test_marketplaces_expose_shared_skills_and_reachable_references() -> None:
    """Both hosts load the same skills and references without external symlinks."""
    manifests = []
    for host, marketplace_path in (
        ("codex", ".agents/plugins/marketplace.json"),
        ("claude", ".claude-plugin/marketplace.json"),
    ):
        marketplace = json.loads((ROOT / marketplace_path).read_text())
        (entry,) = marketplace["plugins"]
        source = entry["source"]
        plugin_root = (ROOT / (source["path"] if host == "codex" else source)).resolve()
        manifest = json.loads((plugin_root / f".{host}-plugin/plugin.json").read_text())
        assert plugin_root == ROOT
        assert marketplace["name"] == entry["name"] == manifest["name"] == "biov"
        assert (plugin_root / manifest["skills"]).resolve() == ROOT / "skills"
        assert (plugin_root / manifest["mcpServers"]).resolve() == ROOT / ".mcp.json"
        manifests.append(manifest)
    assert (
        manifests[0]["version"]
        == manifests[1]["version"]
        == tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    )

    skills = sorted((ROOT / "skills").glob("*/SKILL.md"))
    assert len(skills) == 14
    for skill in skills:
        text = skill.read_text()
        assert re.search(rf"^name: {re.escape(skill.parent.name)}$", text, re.MULTILINE)
        assert re.search(r"^description: .+", text, re.MULTILINE)
        references = re.findall(r"\]\(([^)]+)\)", text)
        references += re.findall(r"`(references/[^`]+)`", text)
        for reference in references:
            if "://" in reference or reference.startswith("#"):
                continue
            target = (skill.parent / reference.split("#", 1)[0]).resolve()
            assert target.is_relative_to(ROOT), (skill, reference)
            assert target.exists(), (skill, reference)


def test_lab_protocols_links_to_publishers_without_bundled_full_text() -> None:
    """Distribute discovery guidance, without redistributing publisher collections."""
    skill_root = ROOT / "skills/lab-protocols"
    assert not any(path.is_file() for path in (skill_root / "references").rglob("*"))
    links = set(
        re.findall(r"\]\((https://[^)]+)\)", (skill_root / "SKILL.md").read_text())
    )
    assert {
        "https://www.addgene.org/protocols/",
        "https://www.thermofisher.com/us/en/home/references/protocols.html",
        "https://www.protocols.io/",
        "https://apidoc.protocols.io/",
    } <= links


def test_declared_plugin_mcp_starts_and_exposes_context() -> None:
    """Launch the manifest command and inspect its actual MCP initialization."""
    config = json.loads((ROOT / ".mcp.json").read_text())["mcpServers"]["biov"]

    async def check() -> None:
        parameters = StdioServerParameters(
            **config,
            env={
                **os.environ,
                "PATH": f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}",
            },
        )
        with anyio.fail_after(30):
            async with stdio_client(parameters) as (reader, writer):
                async with ClientSession(reader, writer) as session:
                    initialized = await session.initialize()
                    assert initialized.instructions is not None
                    assert "scientific-software" in initialized.instructions
                    assert "biov exec" in initialized.instructions
                    assert {t.name for t in (await session.list_tools()).tools} == {
                        "parse_identifiers",
                        "resolve_identifiers",
                        "run_analysis",
                        "inspect_analysis",
                    }
                    result = await session.read_resource("identifiers://uniprot")
                    content = result.contents[0]
                    assert isinstance(content, TextResourceContents)
                    assert json.loads(content.text)["prefix"] == "uniprot"
                    tool_result = await session.call_tool(
                        "resolve_identifiers",
                        {"uri": "identifiers://uniprot"},
                    )
                    embedded = tool_result.content[0]
                    assert embedded.type == "resource"
                    assert embedded.resource == content

    anyio.run(check)
