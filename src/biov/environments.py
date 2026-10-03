"""Install pinned Pixi and scientific environments from native locked manifests."""

import hashlib
import json
import os
import platform
import re
import shutil
import ssl
import subprocess  # noqa: S404 - invoke the pinned package manager with argv
import tarfile
import tempfile
import tomllib
import warnings
import zipfile
from pathlib import Path
from urllib.request import urlopen

import certifi

from .config import settings

_ASSETS = Path(__file__).parent / "assets"
with (_ASSETS / "environments" / "pyproject.toml").open("rb") as _source:
    PIXI_VERSION = tomllib.load(_source)["tool"]["pixi"]["workspace"][
        "requires-pixi"
    ].removeprefix("==")
# Operator override naming an existing Pixi executable to use instead of the
# managed copy; a Settings field of the same name would read this variable.
PIXI_BIN = "BIOV_PIXI_BIN"
# Pixi task setup runs after installation for environments that declare it.
PREPARE_TASK = "prepare"
# A completed preparation is recorded inside the environment it prepared, so the
# record shares that prefix's lifecycle: Pixi cleaning, removing or reinstalling
# the environment invalidates it with the prefix, and a later run can still tell
# "already prepared" from "installed but unusable".
PREPARED_RECORD = ".biov-prepared"
# Pixi resolves platform-scoped tables for the current platform, where a
# ``target.<platform>`` declaration overrides the unscoped one in the same table.
_PLATFORM_TARGETS = {
    ("Linux", "x86_64"): "linux-64",
    ("Linux", "amd64"): "linux-64",
    ("Linux", "aarch64"): "linux-aarch64",
    ("Linux", "arm64"): "linux-aarch64",
    ("Linux", "ppc64le"): "linux-ppc64le",
    ("Linux", "s390x"): "linux-s390x",
    ("Darwin", "x86_64"): "osx-64",
    ("Darwin", "arm64"): "osx-arm64",
    ("Windows", "amd64"): "win-64",
    ("Windows", "x86_64"): "win-64",
    ("Windows", "arm64"): "win-arm64",
}
# SHA-256 digests published with prefix-dev/pixi v0.81.0 release assets.
_RELEASES = {
    ("Darwin", "arm64"): (
        "pixi-aarch64-apple-darwin.tar.gz",
        "f4e32ea91970d4e11739488817979a5f2c6ebbb9cedb0d6dea74b2b790b272dc",
    ),
    ("Darwin", "x86_64"): (
        "pixi-x86_64-apple-darwin.tar.gz",
        "9859588ba57f390b5c77d56b2654fab37bc00e10952da25e3efd6e3434557e5a",
    ),
    ("Linux", "x86_64"): (
        "pixi-x86_64-unknown-linux-musl.tar.gz",
        "7aa3ec39aecceff9062fa2ed4d42cbaa0bdc25ddea727d048e061cf188d434f6",
    ),
    ("Linux", "aarch64"): (
        "pixi-aarch64-unknown-linux-musl.tar.gz",
        "9f8d2113fe9dc01788a65f5c2acec34fa56b1193461a5c3e9a775d6d2d621bcb",
    ),
    ("Windows", "amd64"): (
        "pixi-x86_64-pc-windows-msvc.zip",
        "1fc82219c96d539e6a856bd1bb2643ec9e1753f5b3af61b0c3158921aca89de3",
    ),
    ("Windows", "arm64"): (
        "pixi-aarch64-pc-windows-msvc.zip",
        "4981a25df1a26389712bc1c27992635cff629a635b083755d71938ebe92d5e89",
    ),
}


def managed_pixi_path() -> Path:
    """Return where BioV publishes its own pinned Pixi copy.

    Returns:
        Managed executable path; the file may not exist yet.
    """
    return (
        settings.environment_root.expanduser().resolve()
        / f"pixi-{PIXI_VERSION}"
        / ("pixi.exe" if os.name == "nt" else "pixi")
    )


def _reported_pixi_version(executable: Path) -> str | None:
    """Return the version line a candidate Pixi executable prints.

    Args:
        executable: Candidate Pixi binary.

    Returns:
        Reported ``--version`` output, or None when the candidate cannot run.
    """
    try:
        result = subprocess.run(  # noqa: S603
            [str(executable), "--version"], check=True, capture_output=True, text=True
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout.strip()


def _matching_system_pixi() -> Path | None:
    """Find an already installed Pixi that is exactly the pinned version.

    A Pixi on ``PATH`` with another version, or one that cannot report a
    version, is ignored with a warning so the reason stays visible; the caller
    then falls back to BioV's own pinned copy.

    Returns:
        Matching executable, or None when no suitable Pixi is on ``PATH``.
    """
    found = shutil.which("pixi")
    if found is None:
        return None
    candidate = Path(found).resolve()
    reported = _reported_pixi_version(candidate)
    if reported != f"pixi {PIXI_VERSION}":
        warnings.warn(
            f"Ignoring pixi at {candidate}: expected pixi {PIXI_VERSION}, found "
            f"{reported!r}; BioV falls back to its own pinned copy",
            RuntimeWarning,
            stacklevel=3,
        )
        return None
    return candidate


def pixi_path() -> Path:
    """Return the Pixi executable BioV runs, without installing anything.

    Resolution order is the explicit ``BIOV_PIXI_BIN`` override, a ``PATH``
    Pixi of the pinned version, then BioV's own pinned copy under
    ``environment_root``. The override is a path, or a command name resolved on
    ``PATH``, and is verified by setup like every other candidate.

    Returns:
        Absolute executable path, or the managed path that setup would install.
    """
    override = os.environ.get(PIXI_BIN)
    if override:
        return Path(shutil.which(override) or Path(override).expanduser()).resolve()
    system = _matching_system_pixi()
    if system is not None:
        return system
    return managed_pixi_path()


def pixi_command() -> tuple[str, ...]:
    """Return the argv prefix that invokes the resolved Pixi.

    Returns:
        Command parts for a passthrough alias, to extend with native arguments.
    """
    return (str(pixi_path()),)


def manifest_source() -> Path:
    """Return the configured, packaged, or source-checkout manifest."""
    if settings.environment_manifest is not None:
        return settings.environment_manifest.expanduser().resolve(strict=True)
    return _ASSETS / "environments" / "pyproject.toml"


def _pixi_manifest() -> dict[str, object]:
    """Return the ``tool.pixi`` table of the manifest setup would install.

    Returns:
        Parsed Pixi configuration, or an empty mapping when it declares none.
    """
    with manifest_source().open("rb") as source:
        document = tomllib.load(source)
    section: object = document
    for key in ("tool", "pixi"):
        if not isinstance(section, dict):
            return {}
        section = section.get(key, {})
    return section if isinstance(section, dict) else {}


def declared_environments() -> dict[str, tuple[str, ...]]:
    """Return every Pixi environment the manifest declares, with its features.

    Environments and their scientific dependencies live in the manifest.
    Native feature lists and tables
    with a ``features`` list both select the same features. Setup installs an
    environment only when it appears here.

    Returns:
        Environment name to the features it selects, in declaration order.

    Raises:
        ValueError: If an environment name or feature list is invalid.
    """
    records = _pixi_manifest().get("environments", {})
    if not isinstance(records, dict):
        return {}
    declared: dict[str, tuple[str, ...]] = {}
    for name, record in records.items():
        if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name) is None:
            raise ValueError(f"Invalid Pixi environment name: {name!r}")
        if isinstance(record, list):
            features = record
        elif isinstance(record, dict):
            features = record.get("features", [])
        else:
            features = None
        if not isinstance(features, list) or not all(
            isinstance(feature, str) for feature in features
        ):
            raise ValueError(
                f"[tool.pixi.environments.{name}] must declare a features list"
            )
        declared[name] = tuple(features)
    return declared


def _platform_target() -> str | None:
    """Return the Pixi platform name of this host.

    Returns:
        Platform name such as ``linux-64``, or None when it is unknown so only
        unscoped declarations apply instead of guessing another platform.
    """
    return _PLATFORM_TARGETS.get((platform.system(), platform.machine().lower()))


def _platform_tasks(table: object) -> dict[object, object] | None:
    """Return one manifest table's tasks as this platform resolves them.

    Args:
        table: Workspace or feature table from ``tool.pixi``.

    Returns:
        Task mapping with this platform's ``target`` overrides applied, or None
        when the table declares no task for this platform.
    """
    if not isinstance(table, dict):
        return None
    declared = table.get("tasks")
    tasks: dict[object, object] = dict(declared) if isinstance(declared, dict) else {}
    targets = table.get("target")
    target = _platform_target()
    if target is not None and isinstance(targets, dict):
        record = targets.get(target)
        if isinstance(record, dict) and isinstance(record.get("tasks"), dict):
            tasks.update(record["tasks"])
    return tasks or None


def _prepare_declarations(name: str) -> tuple[object, ...]:
    """Return the ``prepare`` task declarations one environment resolves to.

    Source-only tools ship their native build or checkout steps in Pixi's own
    task table, so preparation is derived from the manifest instead of
    hard-coded environment names. A feature task overrides the workspace
    default, and the default applies only to environments that do not opt out
    of the default feature. Platform-scoped ``target`` declarations override the
    unscoped task in the same table, matching Pixi's own resolution.

    Args:
        name: Pixi environment name, declared or not.

    Returns:
        Task declarations, empty when the environment declares no preparation.
    """
    declared = declared_environments()
    if name not in declared:
        return ()
    manifest = _pixi_manifest()
    features = manifest.get("feature", {})
    declarations: list[object] = []
    for feature in declared[name]:
        record = features.get(feature) if isinstance(features, dict) else None
        tasks = _platform_tasks(record)
        if tasks is not None and PREPARE_TASK in tasks:
            declarations.append(tasks[PREPARE_TASK])
    records = manifest.get("environments", {})
    record = records.get(name) if isinstance(records, dict) else None
    opted_out = isinstance(record, dict) and bool(
        record.get("no-default-feature", False)
    )
    default_tasks = _platform_tasks(manifest)
    if not opted_out and default_tasks is not None and PREPARE_TASK in default_tasks:
        declarations.append(default_tasks[PREPARE_TASK])
    return tuple(declarations)


def prepare_environments() -> tuple[str, ...]:
    """Return the environments that declare a reachable ``prepare`` task.

    Returns:
        Environment names, in manifest declaration order.
    """
    return tuple(
        name for name in declared_environments() if _prepare_declarations(name)
    )


def prepare_command(name: str) -> str | None:
    """Return the canonical ``prepare`` task one environment resolves to.

    Args:
        name: Pixi environment name, declared or not.

    Returns:
        Canonical task text, or None when the environment declares no preparation.
    """
    declarations = _prepare_declarations(name)
    if not declarations:
        return None
    return json.dumps(declarations, sort_keys=True, default=str)


def environment_prefixes(executable: str, manifest: Path) -> dict[str, Path]:
    """Ask Pixi where it installs each declared environment.

    The record of a completed preparation lives inside that environment, so its
    path has to come from Pixi instead of a layout BioV would have to keep
    guessing.

    Args:
        executable: Resolved Pixi executable.
        manifest: Manifest the environments are installed from.

    Returns:
        Environment name to prefix path; empty when Pixi cannot report them, so
        callers treat an unknown prefix as unprepared instead of skipping work.
    """
    try:
        result = subprocess.run(  # noqa: S603
            [
                executable,
                "info",
                "--no-config",
                "--json",
                "--manifest-path",
                str(manifest),
            ],
            capture_output=True,
            check=True,
            text=True,
        )
        records = json.loads(result.stdout)["environments_info"]
    except (OSError, KeyError, TypeError, ValueError, subprocess.SubprocessError):
        # Prefix discovery stays advisory: without it preparation runs again
        # rather than being skipped on a stale record.
        return {}
    prefixes: dict[str, Path] = {}
    for record in records:
        if not isinstance(record, dict):
            continue
        name, prefix = record.get("name"), record.get("prefix")
        if isinstance(name, str) and isinstance(prefix, str) and prefix:
            prefixes[name] = Path(prefix)
    return prefixes


def prepared_marker(prefix: Path) -> Path:
    """Return where a completed preparation of one environment is recorded.

    Args:
        prefix: Environment prefix Pixi installed.

    Returns:
        Marker path inside that environment; it exists only after preparation.
    """
    return prefix / PREPARED_RECORD


def _lock_digest(manifest: Path) -> str:
    """Return the digest of the lock an environment is installed from.

    Args:
        manifest: Manifest the environment is installed from.

    Returns:
        Lock SHA-256, or ``-`` when the manifest ships no lock.
    """
    lock = manifest.with_name("pixi.lock")
    if not lock.is_file():
        return "-"
    with lock.open("rb") as content:
        return hashlib.file_digest(content, "sha256").hexdigest()


def _preparation_record(manifest: Path, name: str) -> str:
    """Return the marker text a completed preparation of one environment writes.

    The record invalidates itself: it changes with the environment name, with
    the ``prepare`` task that environment resolves to, and with the lock. The
    prefix it is written into supplies the remaining lifecycle.

    Args:
        manifest: Manifest the environment is installed from.
        name: Declared Pixi environment name.

    Returns:
        Marker text for the current inputs.
    """
    task = hashlib.sha256((prepare_command(name) or "").encode()).hexdigest()
    return f"environment = {name}\nprepare = {task}\nlock = {_lock_digest(manifest)}\n"


def preparation_recorded(manifest: Path, name: str, prefix: Path | None) -> bool:
    """Report whether one environment was already prepared from these inputs.

    Args:
        manifest: Manifest the environment is installed from.
        name: Declared Pixi environment name.
        prefix: Environment prefix Pixi reported, or None when it reported none.

    Returns:
        True only when the environment still holds a marker matching the current
        task and lock; a cleaned or reinstalled prefix no longer does.
    """
    if prefix is None:
        return False
    try:
        recorded = prepared_marker(prefix).read_text(encoding="utf-8")
    except OSError:
        return False
    return recorded == _preparation_record(manifest, name)


def record_preparation(manifest: Path, name: str, prefix: Path | None) -> None:
    """Record one successful preparation so later runs can skip repeating it.

    Args:
        manifest: Manifest the environment was installed from.
        name: Declared Pixi environment name.
        prefix: Environment prefix Pixi reported. Without one, or when it is
            already gone, nothing is recorded, so the next run prepares again
            instead of skipping work on a record that cannot be reused.
    """
    if prefix is None or not prefix.is_dir():
        return
    prepared_marker(prefix).write_text(
        _preparation_record(manifest, name), encoding="utf-8"
    )


def provision_environment(name: str, *, executable: str | None = None) -> Path:
    """Install a declared locked environment and prepare it when it declares one.

    On-demand provisioning mirrors explicit setup: installation reads the lock
    and aborts when the lock is stale rather than re-solving, and the
    environment's ``prepare`` task runs only once because its completion is
    recorded inside the environment it prepared, so a cleaned or reinstalled
    prefix prepares again while a stale lock or a failed preparation surfaces as
    the manager's own failure. Preparation is recorded after it succeeds, so an
    interrupted run is retried instead of hidden by a marker.

    Args:
        name: Declared Pixi environment to make ready for execution.
        executable: Already resolved Pixi executable; the manager is resolved
            here only when the caller has not resolved it yet.

    Returns:
        Manifest path the environment was installed from.
    """
    manifest = environment_manifest(prepare=True)
    pixi = str(pixi_path()) if executable is None else executable
    common = ["--no-config", "--manifest-path", str(manifest)]
    subprocess.run(  # noqa: S603
        [pixi, "install", *common, "--locked", "--environment", name],
        check=True,
    )
    if name in prepare_environments():
        prefix = environment_prefixes(pixi, manifest).get(name)
        if not preparation_recorded(manifest, name, prefix):
            subprocess.run(  # noqa: S603
                [pixi, "run", *common, "--frozen", "--environment", name, PREPARE_TASK],
                check=True,
            )
            record_preparation(manifest, name, prefix)
    return manifest


def resolve_environment(tool: str) -> str:
    """Require a Pixi environment the manifest declares.

    Args:
        tool: Pixi environment name to install.

    Returns:
        Declared Pixi environment name.

    Raises:
        ValueError: If the argument names no declared environment.
    """
    declared = declared_environments()
    if tool not in declared:
        choices = ", ".join(sorted(declared)) or "none"
        raise ValueError(
            f"Unknown Pixi environment '{tool}'; "
            f"{manifest_source()} declares: {choices}"
        )
    return tool


def environment_manifest(*, prepare: bool = False) -> Path:
    """Locate an explicit project or a writable copy of the bundled manifest.

    The source checkout and wheel share the same pyproject.toml and pixi.lock.
    Bundled workspaces are keyed by content so upgrades do not modify old ones.

    Returns:
        Manifest path; only setup creates the bundled workspace.

    Raises:
        OSError: If the workspace cannot be published.
    """
    if settings.environment_manifest is not None:
        return settings.environment_manifest.expanduser().resolve(strict=True)
    source = manifest_source().parent
    files = {
        name: (source / name).read_bytes() for name in ("pyproject.toml", "pixi.lock")
    }
    digest = hashlib.sha256(b"\0".join(files.values())).hexdigest()
    root = settings.environment_root.expanduser().resolve()
    workspace = root / "workspaces" / digest
    if prepare and not workspace.exists():
        workspace.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=workspace.parent) as staging:
            package = Path(staging) / "workspace"
            package.mkdir()
            for name, content in files.items():
                (package / name).write_bytes(content)
            try:
                os.rename(package, workspace)
            except OSError:
                if not workspace.is_dir():
                    raise
    return workspace / "pyproject.toml"


def _publish_manager(package: Path, executable: Path) -> None:
    """Publish a verified manager package over any stale install directory.

    Args:
        package: Verified staging tree containing the manager binary.
        executable: Managed executable path the package is published as.

    Raises:
        OSError: If the package is not published and no manager is installed.
    """
    parent = executable.parent
    if parent.exists() and not executable.is_file():
        # A partial installation, or an operator who deleted only the binary,
        # leaves a directory that the atomic rename cannot replace; clear it so
        # setup repairs itself instead of failing on a non-empty directory.
        if parent.is_dir():
            shutil.rmtree(parent)
        else:
            parent.unlink()
    try:
        os.rename(package, parent)
    except OSError:
        if not executable.is_file():
            raise


def setup_environment(
    tool: str | None = None,
    *,
    archive: Path | None = None,
    all_environments: bool = False,
    update_lock: bool = False,
) -> Path:
    """Install Pixi, optionally installing one declared locked environment.

    An already installed Pixi of the pinned version is reused: an explicit
    ``BIOV_PIXI_BIN`` override, then a matching ``pixi`` on ``PATH``, then the
    managed copy. Only when none is usable is the pinned release downloaded.
    A ``BIOV_PIXI_BIN`` override that cannot be executed is an error rather
    than a silent fallback, because the operator asked for it explicitly.
    Failing to publish the managed copy surfaces as an ``OSError``.

    Args:
        tool: Pixi environment name, or None to install only
            the manager.
        archive: Official release archive supplied for offline manager setup.
        all_environments: Install every environment declared in the manifest.
        update_lock: Explicitly regenerate an operator-selected project's lock.

    Returns:
        Installed manager path, or the environment's manifest path.

    Raises:
        ValueError: If the environment, execution host, platform or checksum is invalid.
    """
    if settings.execution_host is not None:
        raise ValueError(
            "Run setup on the execution host with a local BioV configuration"
        )
    if tool is not None and all_environments:
        raise ValueError("Choose one environment name or --all")
    if update_lock and settings.environment_manifest is None:
        raise ValueError(
            "Lock updates require BIOV_ENVIRONMENT_MANIFEST pointing to a project"
        )
    environment = resolve_environment(tool) if tool is not None else None
    if (
        (all_environments or environment is not None)
        and settings.environment_manifest is None
        and (platform.system(), platform.machine().lower()) != ("Linux", "x86_64")
    ):
        raise ValueError(
            "Bundled scientific environments require linux-64; select a compatible BIOV_ENVIRONMENT_MANIFEST or use a Linux execution host"
        )
    executable = pixi_path()
    root = settings.environment_root.expanduser().resolve()
    if not executable.is_file():
        override = os.environ.get(PIXI_BIN)
        if override:
            raise ValueError(
                f"{PIXI_BIN} does not name an executable Pixi: {executable}"
            )
        executable = managed_pixi_path()
        release = _RELEASES.get((platform.system(), platform.machine().lower()))
        if release is None:
            raise ValueError("No bundled Pixi release for this platform")
        filename, digest = release
        root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".pixi-", dir=root) as staging:
            source = archive
            if source is None:
                source = Path(staging) / filename
                url = f"https://github.com/prefix-dev/pixi/releases/download/v{PIXI_VERSION}/{filename}"
                context = ssl.create_default_context()
                context.load_verify_locations(certifi.where())
                with (
                    urlopen(url, timeout=60, context=context) as response,
                    source.open("xb") as output,
                ):
                    shutil.copyfileobj(response, output)
            with source.open("rb") as content:
                if hashlib.file_digest(content, "sha256").hexdigest() != digest:
                    raise ValueError("Pixi archive SHA-256 mismatch")
            package = Path(staging) / "package"
            if filename.endswith(".zip"):
                with zipfile.ZipFile(source) as zipped:
                    zipped.extractall(package)  # noqa: S202 - pinned archive digest verified above
            else:
                with tarfile.open(source) as tar:
                    tar.extractall(package, filter="data")
            shutil.copyfile(_ASSETS / "pixi-LICENSE", package / "LICENSE")
            subprocess.run([str(package / executable.name), "--version"], check=True)  # noqa: S603
            _publish_manager(package, executable)
    result = subprocess.run(  # noqa: S603
        [str(executable), "--version"], check=True, capture_output=True, text=True
    )
    if result.stdout.strip() != f"pixi {PIXI_VERSION}":
        raise ValueError(
            f"Expected pixi {PIXI_VERSION}, got {result.stdout.strip()!r} from {executable}"
        )
    if update_lock or all_environments or environment is not None:
        manifest = environment_manifest(prepare=True)
        common = ["--no-config", "--manifest-path", str(manifest)]
        if update_lock:
            subprocess.run([str(executable), "lock", *common], check=True)  # noqa: S603
        if all_environments or environment is not None:
            selection = (
                ["--environment", environment] if environment is not None else ["--all"]
            )
            subprocess.run(  # noqa: S603
                [str(executable), "install", *common, "--locked", *selection],
                check=True,
            )
            installed = (
                tuple(declared_environments()) if all_environments else (environment,)
            )
            preparations = set(prepare_environments())
            prefixes: dict[str, Path] | None = None
            for name in installed:
                if name in preparations:
                    subprocess.run(  # noqa: S603
                        [
                            str(executable),
                            "run",
                            *common,
                            "--as-is",
                            "--environment",
                            name,
                            PREPARE_TASK,
                        ],
                        check=True,
                    )
                    if prefixes is None:
                        # Ask the manager once where this project's records go.
                        prefixes = environment_prefixes(str(executable), manifest)
                    # An explicit setup still prepares unconditionally; the
                    # record only spares the next exec from repeating minutes
                    # of work that already succeeded.
                    record_preparation(manifest, name, prefixes.get(name))
        return manifest
    return executable
