import os
from pathlib import Path

import fsspec.config
from pydantic_settings import BaseSettings


class Settings(BaseSettings, env_prefix="BIOV_"):
    home: Path | None = None
    cache_http: bool = True


settings = Settings()

appname = __package__ or "biov"

if settings.home is not None:
    cache_storage = str(settings.home)
elif (xdg_cache_home := os.getenv("XDG_CACHE_HOME")) is not None:
    cache_storage = os.path.join(xdg_cache_home, appname)
else:
    from platformdirs import PlatformDirs

    cache_storage = PlatformDirs(appname=appname).user_cache_dir

fsspec.config.conf["filecache"] = {"cache_storage": cache_storage}
