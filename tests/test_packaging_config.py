"""Keep the upstream mixed Rust/Python distribution contract explicit."""

import configparser
import fnmatch
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _pyproject() -> dict:
    with (ROOT / "pyproject.toml").open("rb") as stream:
        return tomllib.load(stream)


def test_one_distribution_uses_upstream_rust_bin_and_extension() -> None:
    """One upstream backend builds both mandatory native artifacts."""
    config = _pyproject()
    build = config["build-system"]
    assert build["build-backend"] == "setuptools.build_meta"
    assert "backend-path" not in build
    assert any(
        requirement.startswith("setuptools-rust") for requirement in build["requires"]
    )
    assert any(
        requirement.startswith("setuptools>=") for requirement in build["requires"]
    )
    assert "biov" not in config["project"].get("scripts", {})

    (extension,) = config["tool"]["setuptools-rust"]["ext-modules"]
    (binary,) = config["tool"]["setuptools-rust"]["bins"]
    assert extension["target"] == "biov._native"
    assert extension["path"] == "crates/biov-python/Cargo.toml"
    assert extension.get("binding", "PyO3") == "PyO3"
    assert "extension-module" in extension["features"]
    assert binary["target"] == "biov"
    assert binary["path"] == "crates/biov-cli/Cargo.toml"
    for artifact in (extension, binary):
        assert "--locked" in artifact["cargo-manifest-args"]
        assert artifact.get("optional", False) is False

    with (ROOT / binary["path"]).open("rb") as stream:
        cargo = tomllib.load(stream)
    assert [target["name"] for target in cargo["bin"]] == ["biov"]


def test_existing_abi3_contract_has_matching_wheel_configuration() -> None:
    """Setuptools's stable-ABI wheel tag matches the Rust PyO3 feature."""
    config = configparser.ConfigParser()
    config.read(ROOT / "setup.cfg")
    assert config["bdist_wheel"]["py_limited_api"] == "cp312"
    with (ROOT / "Cargo.toml").open("rb") as stream:
        cargo = tomllib.load(stream)
    assert "abi3-py312" in cargo["workspace"]["dependencies"]["pyo3"]["features"]


def test_python_package_data_covers_every_source_asset() -> None:
    """Registry, environments, type hints and CRISPR resources stay bundled."""
    config = _pyproject()["tool"]["setuptools"]
    assert config["package-dir"][""] == "src"
    find = config["packages"]["find"]
    assert find["where"] == ["src"]
    assert set(find["include"]) == {"biov*", "crisprprimer*"}
    for package in ("biov", "crisprprimer"):
        root = ROOT / "src" / package
        patterns = config["package-data"][package]
        for path in root.rglob("*"):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            if path.suffix in {".py", ".pyc", ".so", ".pyd", ".dylib"}:
                continue
            relative = path.relative_to(root).as_posix()
            assert any(
                fnmatch.fnmatchcase(relative, pattern) for pattern in patterns
            ), relative


def test_sdist_keeps_the_workspace_and_excludes_host_plugin_configuration() -> None:
    """The source distribution can rebuild without copying agent internals."""
    manifest = (ROOT / "MANIFEST.in").read_text()
    assert "include Cargo.toml Cargo.lock rust-toolchain.toml uv.lock" in manifest
    for directory in ("crates", "docs", "tests", "scripts", "src/biov/assets"):
        assert f"graft {directory}" in manifest.splitlines()
    for directory in (
        ".agents",
        ".claude-plugin",
        ".codex-plugin",
        ".github",
        "skills",
    ):
        assert f"prune {directory}" in manifest.splitlines()
    assert "exclude tests/test_plugin_distribution.py" in manifest.splitlines()
    assert "setup.cfg" in manifest


def test_python_distribution_and_native_workspace_versions_match() -> None:
    """The single public command and Python distribution share one release version."""
    with (ROOT / "Cargo.toml").open("rb") as stream:
        cargo = tomllib.load(stream)
    assert (
        _pyproject()["project"]["version"] == cargo["workspace"]["package"]["version"]
    )
