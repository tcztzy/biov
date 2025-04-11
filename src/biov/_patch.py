import os

import fsspec.config
import pandas.io.common

appname = __package__ or "biov"

if (crisprprimer_home := os.getenv("BIOV_HOME")) is not None:
    cache_storage = crisprprimer_home
elif (xdg_cache_home := os.getenv("XDG_CACHE_HOME")) is not None:
    cache_storage = os.path.join(xdg_cache_home, appname)
else:
    from platformdirs import PlatformDirs

    cache_storage = PlatformDirs(appname=appname).user_cache_dir

fsspec.config.conf["filecache"] = {"cache_storage": cache_storage}

# if pandas dependencies upgrade to 3.0, this patch could be removed
if not pandas.io.common.is_fsspec_url("filecache::s3://some-bucket/some-file"):  # type: ignore
    import re

    _FSSPEC_URL_PATTERN = re.compile(
        r"^[A-Za-z][A-Za-z0-9+\-+.]*(::[A-Za-z0-9+\-+.]+)*://"
    )

    def is_fsspec_url(url):
        return (
            isinstance(url, str)
            and bool(_FSSPEC_URL_PATTERN.match(url))
            and not url.startswith(("http://", "https://"))
        )

    pandas.io.common.is_fsspec_url = is_fsspec_url  # type: ignore
