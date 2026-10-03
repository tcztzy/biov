"""Configuration module for BioV."""

import argparse
import os
import sys
from pathlib import Path

import fsspec.config
from platformdirs import user_config_path, user_data_path
from pydantic import AnyHttpUrl, Field, field_validator
from pydantic_settings import (
    BaseSettings,
    InitSettingsSource,
    PydanticBaseSettingsSource,
    TomlConfigSettingsSource,
)

CONFIG_FILE_ENV_VAR = "BIOV_CONFIG"
CONFIG_FILE_OPTION = "--config"


def get_default_config_path() -> Path:
    """Get the default application TOML path.

    Returns:
        Platform-specific user configuration file for BioV.
    """
    return user_config_path("biov") / "config.toml"


class ConfigFileError(ValueError):
    """An explicitly selected application TOML file is unusable."""


def _unusable_config_file(reason: str) -> str:
    """Describe one unusable application TOML selection.

    Returns:
        Message naming both places that select the file.
    """
    return (
        f"{CONFIG_FILE_ENV_VAR} (or the CLI {CONFIG_FILE_OPTION} option) must name "
        f"an existing TOML file: {reason}"
    )


def _validated_config_file(selected: object) -> Path:
    """Resolve an explicitly selected application TOML file.

    Args:
        selected: Value of BIOV_CONFIG or the CLI --config option.

    Returns:
        Absolute path of the selected file.

    Raises:
        ConfigFileError: If the selection is empty, missing, or not a file.
    """
    text = "" if selected is None else str(selected).strip()
    if not text:
        raise ConfigFileError(_unusable_config_file("the value is empty"))
    path = Path(text).expanduser().resolve()
    if not path.exists():
        raise ConfigFileError(_unusable_config_file(f"{path} does not exist"))
    if path.is_dir():
        raise ConfigFileError(_unusable_config_file(f"{path} is a directory"))
    if not path.is_file():
        raise ConfigFileError(_unusable_config_file(f"{path} is not a regular file"))
    return path


def get_default_home() -> Path:
    """Get default BIOV_HOME.

    Returns:
        $XDG_CACHE_HOME/biov if XDG_CACHE_HOME set else determined by platformdirs
    """
    if (xdg_cache_home := os.getenv("XDG_CACHE_HOME")) is not None:
        return Path(os.path.join(xdg_cache_home, "biov"))
    else:
        from platformdirs import PlatformDirs

        return PlatformDirs(appname="biov").user_cache_path


class Settings(BaseSettings, env_prefix="BIOV_"):
    """Settings.

    Attributes:
        config: Application TOML file, selected by BIOV_CONFIG or --config.
        home: Custom cache directory
        cache_http: Cache file from http or not
        max_file_bytes: Optional download/decompressed-file ceiling in bytes.
        environment_root: Persistent directory for Pixi and bundled scientific environments.
        environment_manifest: Project manifest, or None for the bundled manifest.
        execution_host: SSH host alias, or None for local execution.
        execution_cwd: Default working directory on the execution host.
        ssh_config: Optional OpenSSH configuration file on this host.
        analysis_root: Persistent local analysis outputs and execution records.
        analysis_base_url: Existing HTTP storage endpoint serving analysis_root.
    """

    config: Path = Field(default_factory=get_default_config_path, validate_default=True)
    home: Path = Field(default_factory=get_default_home, validate_default=True)
    cache_http: bool = True
    max_file_bytes: int | None = Field(default=None, gt=0)
    environment_root: Path = Field(
        default_factory=lambda: user_data_path("biov") / "environments"
    )
    environment_manifest: Path | None = None
    execution_host: str | None = Field(default=None, min_length=1)
    execution_cwd: Path | None = None
    ssh_config: Path | None = None
    analysis_root: Path = Field(
        default_factory=lambda: user_data_path("biov") / "results"
    )
    analysis_base_url: AnyHttpUrl | None = None

    @field_validator("analysis_base_url")
    @classmethod
    def validate_analysis_base_url(cls, value: AnyHttpUrl | None) -> AnyHttpUrl | None:
        """Require a storage prefix without query, fragment or credentials.

        Returns:
            Validated storage prefix, or None for local-only file access.

        Raises:
            ValueError: If the URL is not a plain HTTP storage prefix.
        """
        if value is not None and (
            value.query or value.fragment or value.username or value.password
        ):
            raise ValueError("analysis_base_url must be a plain HTTP(S) storage prefix")
        return value

    @field_validator("home", mode="after")
    @classmethod
    def normalize_home(cls, value: Path) -> Path:
        """Expand and resolve the cache directory once when settings load.

        Returns:
            Absolute cache directory.
        """
        return value.expanduser().resolve()

    @field_validator("config", mode="after")
    @classmethod
    def normalize_config(cls, value: Path) -> Path:
        """Expand and resolve the selected application TOML path.

        Returns:
            Absolute application TOML path, whether or not it exists.
        """
        return value.expanduser().resolve()

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls: type[BaseSettings],
        init_settings: PydanticBaseSettingsSource,
        env_settings: PydanticBaseSettingsSource,
        dotenv_settings: PydanticBaseSettingsSource,
        file_secret_settings: PydanticBaseSettingsSource,
    ) -> tuple[PydanticBaseSettingsSource, ...]:
        """Load operator settings without exposing deployment choices to callers.

        An explicitly selected application TOML must be usable; only the
        platform default may be absent. An unusable selection raises
        ``ConfigFileError`` through ``_validated_config_file``.

        Returns:
            Explicit arguments, environment variables, then the application TOML.
        """
        selected = (
            init_settings.init_kwargs.get("config")
            if isinstance(init_settings, InitSettingsSource)
            else None
        )
        if selected is None:
            selected = os.getenv(CONFIG_FILE_ENV_VAR)
        path = (
            get_default_config_path()
            if selected is None
            else _validated_config_file(selected)
        )
        return init_settings, env_settings, TomlConfigSettingsSource(settings_cls, path)


def _export_selection(selected: Path) -> None:
    """Make one selected application TOML the selection child processes inherit.

    An executed analysis script is a separate process that imports BioV again.
    Without this it would read the platform default, or a stale inherited
    ``BIOV_CONFIG``, instead of the file this invocation selected.

    Args:
        selected: Resolved application TOML path this invocation selected.
    """
    os.environ[CONFIG_FILE_ENV_VAR] = str(selected)


# The console entry point imports this module before Typer can select --config.
# Stop at the subcommand so native programs keep their own --config arguments.
if Path(sys.argv[0]).name in {"biov", "biov.exe"}:
    _config_parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    _config_parser.add_argument(CONFIG_FILE_OPTION, type=Path)
    _config_parser.add_argument("command", nargs=argparse.REMAINDER)
    _cli_config, _ = _config_parser.parse_known_args(sys.argv[1:])
    try:
        settings = (
            Settings(config=_cli_config.config)
            if _cli_config.config is not None
            else Settings()
        )
    except ConfigFileError as error:
        sys.stderr.write(f"{error}\n")
        raise SystemExit(2) from error
    if _cli_config.config is not None:
        # The explicit selection outranks an inherited BIOV_CONFIG for this
        # whole invocation, including the processes it starts.
        _export_selection(settings.config)
else:
    settings = Settings()


_installed_cache_home: str | None = None


def _install_fsspec_cache_default() -> None:
    """Use BIOV_HOME as fsspec's file cache unless the operator configured one."""
    global _installed_cache_home
    filecache = fsspec.config.conf.setdefault("filecache", {})
    current = filecache.get("cache_storage")
    if current is not None and current != _installed_cache_home:
        return
    _installed_cache_home = str(settings.home)
    filecache["cache_storage"] = _installed_cache_home


_install_fsspec_cache_default()


def select_config_file(path: Path | str) -> Settings:
    """Re-read the process-wide settings from an explicitly selected TOML file.

    The module-level ``settings`` object keeps its identity while its fields are
    replaced, so modules that imported it observe the selected application
    settings. An unusable file raises ``ConfigFileError`` before anything is
    replaced. The selection is also exported as ``BIOV_CONFIG`` for this
    invocation, so child processes such as an executed analysis script read the
    same application settings.

    Args:
        path: Application TOML file selected by the CLI --config option.

    Returns:
        The reloaded process-wide settings.
    """
    reloaded = Settings(config=Path(path))
    vars(settings).update(vars(reloaded))
    # Replacing the same object's fields is the point of selecting a file after
    # import, and pydantic's stubs declare this field read-only.
    setattr(  # noqa: B010 - no safer alternative for a stub-declared read-only field.
        settings, "__pydantic_fields_set__", set(reloaded.__pydantic_fields_set__)
    )
    _export_selection(settings.config)
    _install_fsspec_cache_default()
    return settings
