from typing import Any

from ..config import settings


def preprocessing(filepath_or_buffer: Any, **kwargs) -> tuple[Any, Any]:
    if isinstance(filepath_or_buffer, str):
        *protocols, path = filepath_or_buffer.split("::")
        # https://github.com/pandas-dev/pandas/pull/60100
        if any([protocol.startswith("tar://") for protocol in protocols]):
            kwargs["compression"] = None
        if (
            settings.cache_http
            and path.startswith(("https://", "http://"))
            and (len(protocols) == 0 or "filecache" != protocols[-1])
        ):
            filepath_or_buffer = "::".join([*protocols, "filecache", path])
    return filepath_or_buffer, kwargs
