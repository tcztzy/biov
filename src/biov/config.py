"""Configuration module for BioV."""

import os
from pathlib import Path

import fsspec.config
from pydantic import Field
from pydantic_settings import BaseSettings

appname = __package__ or "biov"


def get_default_home() -> Path:
    """Get default BIOV_HOME.

    Returns:
        $XDG_CACHE_HOME/biov if XDG_CACHE_HOME set else determined by platformdirs
    """
    if (xdg_cache_home := os.getenv("XDG_CACHE_HOME")) is not None:
        return Path(os.path.join(xdg_cache_home, appname))
    else:
        from platformdirs import PlatformDirs

        return PlatformDirs(appname=appname).user_cache_path


class Settings(BaseSettings, env_prefix="BIOV_"):
    """Settings.

    Attributes:
        home: Custom cache directory
        cache_http: Cache file from http or not
    """

    home: Path = Field(default_factory=get_default_home)
    cache_http: bool = True


settings = Settings()


fsspec.config.conf["filecache"] = {"cache_storage": str(settings.home)}
