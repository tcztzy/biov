#!/usr/bin/env python3
"""Upstream CLI surface fixture, independent of private backend receipts."""

# ruff: noqa: T201, S101 - fixture deliberately exercises native stdio and assertions
import json
import os
import shutil
import sys
from pathlib import Path

home = Path(__file__).parent
args = sys.argv[1:]
uv = "uv" in Path(sys.argv[0]).name
with (home / "global-argv.jsonl").open("a") as out:
    out.write(
        json.dumps(
            {
                "uv": uv,
                "args": args,
                "env": {
                    key: os.environ.get(key)
                    for key in (
                        "PIXI_HOME",
                        "PIXI_CACHE_DIR",
                        "UV_TOOL_DIR",
                        "UV_TOOL_BIN_DIR",
                        "UV_CACHE_DIR",
                        "UV_PYTHON_INSTALL_DIR",
                    )
                },
            }
        )
        + "\n"
    )
if args == ["--version"]:
    print("uv 0.12.19 (fixture)" if uv else "pixi 0.81.0")
    sys.exit()
if (home / "global-fail").exists():
    sys.exit(23)
if uv:
    root = Path(os.environ["UV_TOOL_DIR"])
    bin_dir = Path(os.environ["UV_TOOL_BIN_DIR"])
    selected = root / "goatools"
    assert args[:1] == ["tool"]
    action = args[1]
    if action == "list":
        if selected.is_dir():
            print("goatools v1.6.5 [with statsmodels==0.14.6]\n- goatools")
        sys.exit()
    commands = (
        "goatools",
        "find_enrichment.py",
        "go_plot.py",
        "map_to_slim.py",
        "ncbi_gene_results_to_python.py",
        "plot_go_term.py",
        "wr_hier.py",
    )
    if action == "install":
        assert args[2] == "goatools==1.6.5"
        assert args[args.index("--with") + 1] == "statsmodels==0.14.6"
        assert "--no-config" in args and "--no-python-downloads" in args
        (selected / "bin").mkdir(parents=True, exist_ok=True)
        bin_dir.mkdir(parents=True, exist_ok=True)
        for name in commands:
            target = selected / "bin" / name
            target.write_text("#!/bin/sh\nexit 0\n")
            target.chmod(0o755)
            if not (bin_dir / name).exists():
                (bin_dir / name).symlink_to(target)
        sys.exit()
    if action == "uninstall":
        for name in commands:
            (bin_dir / name).unlink(missing_ok=True)
        shutil.rmtree(selected)
        sys.exit()
else:
    root = Path(os.environ["PIXI_HOME"])
    assert (root / "manifests/pixi-global.toml").is_file()
    assert args[:1] == ["global"]
    action = args[1]
    selected = root / "envs/samtools"
    if action == "list":
        if "--json" in args:
            print(
                json.dumps(
                    [
                        {
                            "name": "samtools",
                            "dependencies": [{"name": "samtools", "version": "1.24"}],
                            "exposed": [
                                {"exposed_name": "samtools", "executable": "samtools"}
                            ],
                        }
                    ]
                    if selected.is_dir()
                    else []
                )
            )
        else:
            print(
                "samtools 1.24\n└── exposes: samtools"
                if selected.is_dir()
                else "No global environments found."
            )
        sys.exit()
    if action == "install":
        assert "samtools==1.24" in args and "--no-shortcuts" in args
        assert args[args.index("--environment") + 1] == "samtools"
        assert args[args.index("--expose") + 1] == "samtools"
        selected.mkdir(parents=True, exist_ok=True)
        native_bin = root / "bin/trampoline_configuration"
        native_bin.mkdir(parents=True, exist_ok=True)
        (native_bin / "trampoline_bin").write_bytes(b"native trampoline")
        target = root / "bin/samtools"
        target.write_bytes(b"native trampoline")
        target.chmod(0o755)
        sys.exit()
    if action == "uninstall":
        (root / "bin/samtools").unlink(missing_ok=True)
        shutil.rmtree(selected)
        sys.exit()
sys.exit(99)
