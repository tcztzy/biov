"""Read identifier-backed BioV artifacts through fsspec."""

from typing import Any

from fsspec import AbstractFileSystem
from fsspec.implementations.local import LocalFileSystem

from . import artifacts
from .registry import DATA_RESOURCE_NAMESPACE_PREFIXES


class BioVFileSystem(AbstractFileSystem):
    """Read declared biological files through their full identifier URIs.

    Downloads, validation, and caching use the same providers as ``biov.path``.
    File metadata lookup can download a file on a cache miss. Directory listing
    and writes are unsupported, and a provider failure is never reported as a
    missing file.
    """

    protocol = tuple(sorted(DATA_RESOURCE_NAMESPACE_PREFIXES))
    local_file = True

    def __init__(self, *, artifact: str | None = None, **kwargs: Any) -> None:
        """Select the artifact representation for this filesystem.

        Args:
            artifact: File kind, or None for the identifier namespace's default.
            **kwargs: Standard fsspec filesystem options.
        """
        super().__init__(**kwargs)
        self.artifact = artifact

    @classmethod
    def _strip_protocol(cls, path: str) -> str:
        """Keep the namespace in the URI for BioV identifier validation.

        Args:
            path: Full identifier URI.

        Returns:
            Unchanged URI, including its scheme.
        """
        return path

    def _open(
        self,
        path: str,
        mode: str = "rb",
        block_size: int | None = None,
        autocommit: bool = True,
        cache_options: dict[str, Any] | None = None,
        **kwargs: Any,
    ) -> Any:
        """Open the original cached file with fsspec's local file implementation.

        Args:
            path: Full identifier URI.
            mode: Binary read mode; fsspec handles text wrapping.
            block_size: Unused; reads use the local file's buffering.
            autocommit: Unused for read-only files.
            cache_options: Unused; BioV owns the artifact cache.
            **kwargs: Standard fsspec open options; no extra buffering is needed.

        Returns:
            Seekable local file object with its real name and byte size.

        Raises:
            ValueError: If a write or update mode is requested.
        """
        if mode != "rb":
            raise ValueError("BioV filesystems support read-only modes")
        artifact = artifacts.path(path, artifact=self.artifact)
        return LocalFileSystem().open(str(artifact.path), mode="rb")

    def info(self, path: str, **kwargs: Any) -> dict[str, Any]:
        """Return metadata for the original provider file.

        Args:
            path: Full identifier URI.
            **kwargs: Standard fsspec metadata options.

        Returns:
            Identifier URI, byte size, and file type for the selected artifact.
        """
        artifact = artifacts.path(path, artifact=self.artifact)
        return {"name": path, "size": artifact.size, "type": "file"}

    def exists(self, path: str, **kwargs: Any) -> bool:
        """Report resolution without masking a provider failure as "not found".

        fsspec's default implementation swallows every exception, which hides
        provider outages and interrupts behind a missing-file result.

        Args:
            path: Full identifier URI.
            **kwargs: Standard fsspec metadata options.

        Returns:
            True for a resolved artifact; False only when the provider reports
            that the identifier has no such file.
        """
        try:
            self.info(path, **kwargs)
        except FileNotFoundError:
            return False
        return True

    def mkdir(self, path: str, create_parents: bool = True, **kwargs: Any) -> None:
        """Reject directory creation, which fsspec would otherwise ignore.

        Raises:
            NotImplementedError: BioV filesystems are read-only.
        """
        raise NotImplementedError(
            "BioV filesystems are read-only: mkdir is unsupported"
        )

    def makedirs(self, path: str, exist_ok: bool = False) -> None:
        """Reject recursive directory creation instead of silently succeeding.

        Raises:
            NotImplementedError: BioV filesystems are read-only.
        """
        raise NotImplementedError(
            "BioV filesystems are read-only: makedirs is unsupported"
        )

    def rmdir(self, path: str) -> None:
        """Reject directory removal instead of silently succeeding.

        Raises:
            NotImplementedError: BioV filesystems are read-only.
        """
        raise NotImplementedError(
            "BioV filesystems are read-only: rmdir is unsupported"
        )
