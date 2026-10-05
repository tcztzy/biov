#!/usr/bin/env python3
"""Protocol fixture; no Hub client, network access or model execution."""

# ruff: noqa: T201, S101 - intentional native CLI protocol stdio and assertions
import json
import os
import pathlib
import sys

home = pathlib.Path(os.environ["BIOV_FAKE_CLIENT_HOME"])
words = sys.argv[1:]
kind = "uv" if words and words[0] == "tool" else "hf"
local_dir = words[words.index("--local-dir") + 1] if "--local-dir" in words else None
with (home / "calls.jsonl").open("a") as stream:
    stream.write(
        json.dumps(
            {
                "backend": kind,
                "args": words,
                "cwd": str(pathlib.Path.cwd()),
                "local_dir": local_dir,
            }
        )
        + "\n"
    )
if kind == "uv":
    assert words[:9] == [
        "tool",
        "run",
        "--no-config",
        "--no-python-downloads",
        "--from",
        "huggingface-hub==2.1.1",
        "--with",
        "httpx2[socks]",
        "hf",
    ]
    words = words[9:]
if words == ["version"]:
    version_file = home / ("uv-version" if kind == "uv" else "version")
    version = version_file.read_text() if version_file.exists() else "2.1.1"
    print(f"✓ hf version\n  version: {version}")
    sys.exit(0)
if words == ["download", "--help"]:
    print("hf download REPO FILE... --repo-type --revision --local-dir")
    sys.exit(0)
assert words[0] == "download"
assert words[words.index("--repo-type") + 1] == "model"
directory = pathlib.Path(words[words.index("--local-dir") + 1])
assert directory.is_absolute()
assert len(words[words.index("--revision") + 1]) == 40
assert not any(directory.iterdir()), "hf download must use fresh empty staging"
if (home / "relative-environment").exists():
    assert pathlib.Path.cwd() == home
    for name in [
        "HF_HOME",
        "HF_TOKEN_PATH",
        "HF_HUB_CACHE",
        "UV_TOOL_DIR",
        "UV_CACHE_DIR",
    ]:
        assert pathlib.Path(os.environ[name]).read_text() == "existing caller state"
files = words[2 : words.index("--repo-type")]
for name in files:
    if (home / "omit-file").exists() and name == (home / "omit-file").read_text():
        continue
    target = directory / name
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(("native model bytes for " + name).encode())
    if (home / "fail-download").exists():
        print("native hf failure", file=sys.stderr)
        sys.exit(23)
if (home / "concurrent-destination").exists():
    destination = pathlib.Path((home / "concurrent-destination").read_text())
    destination.mkdir()
    (destination / "keep.txt").write_bytes(b"concurrent external data")
print("upstream download output")
print("upstream progress", file=sys.stderr)
