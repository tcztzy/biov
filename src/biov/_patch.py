import pandas.io.common

# if pandas dependencies upgrade to 3.0, this patch could be removed
if not pandas.io.common.is_fsspec_url("filecache::s3://some-bucket/some-file"):  # type: ignore
    import re

    # https://github.com/pandas-dev/pandas/pull/61041
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
