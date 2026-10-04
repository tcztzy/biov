"""Check native source distribution boundaries and required build inputs."""


# This build-verification script deliberately asserts distribution invariants.
# ruff: noqa: S101

import re
import sys
import tarfile
from pathlib import Path, PurePosixPath


def main(filename: str) -> None:
    """Reject missing sources, build artifacts and workstation-specific paths."""
    with tarfile.open(filename) as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        paths = {member.name.split("/", 1)[1] for member in members}
        required = {
            "Cargo.toml",
            "Cargo.lock",
            "uv.lock",
            "rust-toolchain.toml",
            "pyproject.toml",
            "MANIFEST.in",
            "setup.cfg",
            "LICENSE",
            "README.md",
            "SPEC.md",
            "src/biov/_native.pyi",
            "src/biov/py.typed",
            "crates/biov-core/src/sequence.rs",
            "crates/biov-python/src/lib.rs",
            "crates/biov-cli/src/main.rs",
            "crates/biov-cli/src/python_bridge.rs",
            "crates/biov-tools/src/installed.rs",
            "src/biov/_bridge.py",
            "crates/biov-core/src/bin/biov-core.rs",
            "tests/test_native_sequence.py",
            "tests/test_native_metrics.py",
            "docs/guides/sequence-contract.md",
        }
        assert required <= paths, required - paths
        forbidden = {"target", ".venv", "__pycache__", ".git", ".pixi", "site", "dist"}
        selected = {
            "AGENTS.md",
            "Cargo.toml",
            "Cargo.lock",
            "uv.lock",
            "LICENSE",
            "README.md",
            "SPEC.md",
            "mkdocs.yml",
            "pyproject.toml",
            "MANIFEST.in",
            "setup.cfg",
            "rust-toolchain.toml",
        }
        for directory in (
            "src/biov",
            "src/crisprprimer",
            "crates",
            "docs",
            "tests",
            "scripts",
        ):
            selected.update(
                path.as_posix()
                for path in Path(directory).rglob("*")
                if path.is_file()
                and not forbidden.intersection(path.parts)
                and path.suffix not in {".so", ".pyc", ".whl"}
            )
        selected.discard("tests/test_plugin_distribution.py")
        # setuptools generates exactly these metadata files; no arbitrary
        # egg-info subtree or other source omissions are silently accepted.
        generated = {"PKG-INFO"} | {
            f"src/biov.egg-info/{name}"
            for name in (
                "PKG-INFO",
                "SOURCES.txt",
                "dependency_links.txt",
                "entry_points.txt",
                "requires.txt",
                "top_level.txt",
            )
        }
        assert paths - generated == selected, (
            f"sdist mismatch: missing={sorted(selected - paths)}, "
            f"unexpected={sorted(paths - selected - generated)}"
        )
        workstation = re.compile(rb"/(?:Users|Volumes)/")
        for member in members:
            path = PurePosixPath(member.name)
            assert not forbidden.intersection(path.parts), path
            assert path.suffix not in {".so", ".pyc", ".whl"}, path
            stream = archive.extractfile(member)
            assert stream is not None
            assert not workstation.search(stream.read()), path


if __name__ == "__main__":
    main(sys.argv[1])
