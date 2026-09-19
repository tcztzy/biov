"""Identifiers.org Compact Identifier parsing and resolution."""

import json
import re
from dataclasses import dataclass
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote, unquote, urlsplit
from urllib.request import Request, urlopen

from anyio.to_thread import run_sync

from .registry import (
    DATA_RESOURCE_NAMESPACE_PREFIXES,
    URI_ACCESSION_SAFE,
    build_namespace_resource_uri,
    namespaces_by_prefix,
)

RESOLVER_BASE_URL = "https://resolver.api.identifiers.org"
IDENTIFIERS_ORG_BASE_URL = "https://identifiers.org"
DEFAULT_RESOLVER_TIMEOUT_SECONDS = 10.0
MAX_PROMPT_CANDIDATES = 50
BARE_IDENTIFIER_NAMESPACE_ALLOWLIST = frozenset({"refseq.gcf"})

_PREFIX = r"(?:(?:[A-Za-z0-9][A-Za-z0-9._-]*)/)?[A-Za-z0-9][A-Za-z0-9._-]*"
_COMPACT_IDENTIFIER = re.compile(rf"^{_PREFIX}:(?!//)\S+$")
_COMPACT_IDENTIFIER_IN_TEXT = re.compile(
    rf"(?<![A-Za-z0-9._/-])(?P<identifier>{_PREFIX}:(?!//)[^\s<>\[\]{{}}\"'`]+)"
)
_IDENTIFIERS_ORG_URL = re.compile(
    r"https?://(?:www\.)?identifiers\.org/[^\s<>\[\]{}\"'`]+",
    flags=re.IGNORECASE,
)
_DATA_RESOURCE_URI_IN_TEXT = re.compile(
    r"(?<![A-Za-z0-9+.-])"
    r"(?P<scheme>[A-Za-z][A-Za-z0-9+.-]*)://"
    r"(?P<accession>[^\s<>\[\]{}\"'`]+)"
)
_IDENTIFIERS_RESOURCE_URI_IN_TEXT = re.compile(
    r"(?<![A-Za-z0-9+.-])identifiers://"
    r"(?P<registry>[A-Za-z0-9][A-Za-z0-9._~-]*):"
    r"(?P<accession>[^\s<>\[\]{}\"'`]+)",
    flags=re.IGNORECASE,
)
_BARE_IDENTIFIER_IN_TEXT = re.compile(
    r"(?<![A-Za-z0-9._/:\-])"
    r"(?P<identifier>[A-Za-z0-9][A-Za-z0-9._-]*)"
    r"(?![A-Za-z0-9._/\-])"
)
_TRAILING_PROSE = ".,;!?。，；！？、"
_UNBALANCED_CLOSERS = {")": "(", "]": "[", "}": "{"}

Resolution = dict[str, Any]


@dataclass(frozen=True, slots=True)
class IdentifierRef:
    """One locally validated identifier with canonical registry coordinates."""

    namespace: str
    accession: str
    compact_id: str
    resource_uri: str
    provider_code: str | None = None


class IdentifierSyntaxError(ValueError):
    """A value is not one exact identifier form supported by the registry."""

    def __init__(self, value: str, detail: str) -> None:
        """Initialize a stable local parsing failure."""
        self.value = value
        self.detail = detail
        super().__init__(f"Invalid identifier {value!r}: {detail}")


class IdentifierResolutionError(RuntimeError):
    """Base class for expected identifiers.org resolution failures."""


class IdentifierNotFoundError(IdentifierResolutionError):
    """An input is not a resolver-valid Compact Identifier."""

    def __init__(self, compact_id: str, detail: str) -> None:
        """Initialize a validation failure.

        Args:
            compact_id: Compact Identifier that failed validation.
            detail: Human-readable reason returned by the resolver.
        """
        self.compact_id = compact_id
        self.detail = detail
        super().__init__(f"Cannot resolve {compact_id!r}: {detail}")


class IdentifierServiceError(IdentifierResolutionError):
    """The identifiers.org service failed or returned an invalid response."""


def _validate_compact_identifier(compact_id: str) -> str:
    """Validate the local Compact Identifier envelope.

    Args:
        compact_id: Candidate Compact Identifier.

    Returns:
        The unchanged identifier.

    Raises:
        IdentifierNotFoundError: If the value cannot be a Compact Identifier.
    """
    if not _COMPACT_IDENTIFIER.fullmatch(compact_id):
        raise IdentifierNotFoundError(
            compact_id, "expected [provider/]namespace:accession"
        )
    return compact_id


def build_identifiers_org_url(compact_id: str) -> str:
    """Build a persistent identifiers.org URL.

    Args:
        compact_id: Explicit Compact Identifier.

    Returns:
        The persistent identifiers.org URL.
    """
    _validate_compact_identifier(compact_id)
    encoded = quote(compact_id, safe=URI_ACCESSION_SAFE)
    return f"{IDENTIFIERS_ORG_BASE_URL}/{encoded}"


def _trim_trailing_prose(candidate: str) -> str:
    """Remove delimiters introduced by surrounding prose.

    Args:
        candidate: Candidate token including possible prose punctuation.

    Returns:
        The token without unbalanced trailing delimiters.
    """
    candidate = candidate.rstrip(_TRAILING_PROSE)
    while candidate and candidate[-1] in _UNBALANCED_CLOSERS:
        closer = candidate[-1]
        opener = _UNBALANCED_CLOSERS[closer]
        if candidate.count(closer) <= candidate.count(opener):
            break
        candidate = candidate[:-1].rstrip(_TRAILING_PROSE)
    return candidate


def _identifier_from_url(url: str) -> str | None:
    """Convert a current or legacy identifiers.org URL to a Compact Identifier.

    Args:
        url: URL captured from prompt text.

    Returns:
        A Compact Identifier, or ``None`` when the URL has no entity path.
    """
    trimmed = _trim_trailing_prose(url)
    path = unquote(urlsplit(trimmed).path).lstrip("/")
    if not path:
        return None
    if ":" not in path:
        namespace, separator, accession = path.partition("/")
        if not separator or not accession:
            return None
        path = f"{namespace}:{accession}"
    path = _trim_trailing_prose(path)
    return path if _COMPACT_IDENTIFIER.fullmatch(path) else None


def parse_identifier(value: str) -> IdentifierRef:
    """Parse one exact identifier without resolver or network access.

    Accepted values are an explicit Compact Identifier, an identifiers.org URL,
    a BioV identifier resource URI, or an accession belonging to the
    curated bare-ID allowlist. This API intentionally rejects prompt prose and
    multiple identifiers; use the MCP ``parse_identifiers`` tool for prompts.

    Args:
        value: Exact identifier reference to parse.

    Returns:
        Canonical registry namespace, accession, Compact ID, and resource URI.

    Raises:
        IdentifierSyntaxError: If the complete value is unsupported or invalid.
    """
    if not isinstance(value, str):
        raise IdentifierSyntaxError(str(value), "expected a string")
    candidate = value.strip()
    if not candidate:
        raise IdentifierSyntaxError(value, "value is empty")

    namespaces = namespaces_by_prefix()

    if _IDENTIFIERS_ORG_URL.fullmatch(candidate):
        compact_id = _identifier_from_url(candidate)
        if compact_id is None:
            raise IdentifierSyntaxError(value, "identifiers.org URL has no valid ID")
        candidate = compact_id

    provider_code = None
    is_compact = False
    if resource_match := _IDENTIFIERS_RESOURCE_URI_IN_TEXT.fullmatch(candidate):
        namespace_token = unquote(resource_match.group("registry")).casefold()
        accession = unquote(resource_match.group("accession"))
    elif resource_match := _DATA_RESOURCE_URI_IN_TEXT.fullmatch(candidate):
        namespace_token = resource_match.group("scheme").casefold()
        if namespace_token not in DATA_RESOURCE_NAMESPACE_PREFIXES:
            raise IdentifierSyntaxError(
                value,
                f"unknown data URI scheme {resource_match.group('scheme')!r}",
            )
        accession = unquote(resource_match.group("accession"))
    elif _COMPACT_IDENTIFIER.fullmatch(candidate):
        prefix_part, _, accession = candidate.partition(":")
        provider_part, separator, namespace_token = prefix_part.rpartition("/")
        provider_code = provider_part if separator else None
        is_compact = True
    else:
        bare_matches = [
            prefix
            for prefix in BARE_IDENTIFIER_NAMESPACE_ALLOWLIST
            if re.fullmatch(namespaces[prefix]["pattern"], candidate) is not None
        ]
        if len(bare_matches) > 1:  # pragma: no cover - allowlist is unambiguous
            raise IdentifierSyntaxError(value, "bare identifier is ambiguous")
        if not bare_matches:
            raise IdentifierSyntaxError(
                value,
                "expected one explicit Compact ID, BioV resource URI, "
                "identifiers.org URL, or allowlisted bare ID",
            )
        namespace_token = bare_matches[0]
        accession = candidate

    namespace = namespaces.get(namespace_token.casefold())
    if namespace is None:
        raise IdentifierSyntaxError(
            value, f"unknown registry namespace {namespace_token!r}"
        )
    if is_compact and namespace["namespaceEmbeddedInLui"]:
        accession = f"{namespace_token}:{accession}"
    if re.fullmatch(namespace["pattern"], accession) is None:
        raise IdentifierSyntaxError(
            value,
            f"accession does not match registry rule for {namespace['prefix']!r}",
        )
    compact_id = (
        accession
        if namespace["namespaceEmbeddedInLui"]
        else f"{namespace['prefix']}:{accession}"
    )
    if provider_code is not None:
        compact_id = f"{provider_code}/{compact_id}"
    return IdentifierRef(
        namespace=namespace["prefix"],
        accession=accession,
        compact_id=compact_id,
        resource_uri=build_namespace_resource_uri(namespace, accession),
        provider_code=provider_code,
    )


def extract_identifier_candidates(prompt: str) -> list[str]:
    """Extract supported identifier reference forms from text.

    Candidates are syntactic only; the MCP tool validates each one against the
    official resolver before returning a resource link. BioV resource URI forms
    and bare accessions use the same local parser as ``parse_identifier``.

    Args:
        prompt: User or model text to inspect.

    Returns:
        First-occurrence ordered, deduplicated candidates, capped at
        ``MAX_PROMPT_CANDIDATES``.
    """
    positioned: list[tuple[int, str]] = []
    for match in _IDENTIFIERS_ORG_URL.finditer(prompt):
        if identifier := _identifier_from_url(match.group()):
            positioned.append((match.start(), identifier))

    for match in _COMPACT_IDENTIFIER_IN_TEXT.finditer(prompt):
        candidate = _trim_trailing_prose(match.group("identifier"))
        if _COMPACT_IDENTIFIER.fullmatch(candidate):
            positioned.append((match.start(), candidate))

    for pattern in (_DATA_RESOURCE_URI_IN_TEXT, _BARE_IDENTIFIER_IN_TEXT):
        for match in pattern.finditer(prompt):
            if pattern is _DATA_RESOURCE_URI_IN_TEXT and match.group(
                "scheme"
            ).casefold() not in DATA_RESOURCE_NAMESPACE_PREFIXES | {"identifiers"}:
                continue
            try:
                reference = parse_identifier(_trim_trailing_prose(match.group()))
            except IdentifierSyntaxError:
                continue
            positioned.append((match.start(), reference.compact_id))

    return list(
        dict.fromkeys(
            candidate for _, candidate in sorted(positioned, key=lambda item: item[0])
        )
    )[:MAX_PROMPT_CANDIDATES]


def _error_detail(body: bytes, fallback: str) -> str:
    """Extract an identifiers.org error message without leaking raw responses.

    Args:
        body: Response body returned by the resolver.
        fallback: Message to use when the body has no structured error.

    Returns:
        A concise upstream error message.
    """
    try:
        response = json.loads(body)
    except (UnicodeDecodeError, json.JSONDecodeError):
        return fallback
    message = response.get("errorMessage") if isinstance(response, dict) else None
    return message if isinstance(message, str) and message else fallback


async def resolve_identifier(
    compact_id: str,
    *,
    timeout_seconds: float = DEFAULT_RESOLVER_TIMEOUT_SECONDS,
) -> Resolution:
    """Resolve an identifier on a worker thread without blocking MCP.

    Args:
        compact_id: Locally validated Compact Identifier.
        timeout_seconds: Per-request timeout.

    Returns:
        The original resolver JSON object.

    Raises:
        ValueError: If the timeout is not positive.
    """
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")
    _validate_compact_identifier(compact_id)
    return await run_sync(_resolve_identifier, compact_id, timeout_seconds)


def _resolve_identifier(compact_id: str, timeout_seconds: float) -> Resolution:
    """Perform one blocking resolver request on a worker thread.

    Args:
        compact_id: Locally validated Compact Identifier.
        timeout_seconds: Per-request timeout.

    Returns:
        The official resolver JSON object.

    Raises:
        IdentifierNotFoundError: If the resolver rejects the identifier.
        IdentifierServiceError: If transport or response decoding fails.
    """
    encoded = quote(compact_id, safe="/:")
    request = Request(  # noqa: S310 - origin is a fixed HTTPS constant
        f"{RESOLVER_BASE_URL}/{encoded}",
        headers={
            "Accept": "application/json",
            "User-Agent": "BioV identifiers.org MCP",
        },
    )
    try:
        with urlopen(request, timeout=timeout_seconds) as response:  # noqa: S310
            body = response.read()
    except HTTPError as error:
        detail = _error_detail(error.read(), f"resolver returned HTTP {error.code}")
        if error.code in {400, 404}:
            raise IdentifierNotFoundError(compact_id, detail) from error
        raise IdentifierServiceError(
            f"identifiers.org resolver returned HTTP {error.code}"
        ) from error
    except (TimeoutError, URLError, OSError) as error:
        raise IdentifierServiceError(
            "identifiers.org resolver is unavailable"
        ) from error

    try:
        resolution = json.loads(body)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise IdentifierServiceError(
            "identifiers.org resolver returned invalid JSON"
        ) from error
    if not isinstance(resolution, dict):
        raise IdentifierServiceError("identifiers.org resolver returned invalid JSON")
    if error_message := resolution.get("errorMessage"):
        detail = str(error_message)
        raise IdentifierNotFoundError(compact_id, detail)
    if not isinstance(resolution.get("payload"), dict):
        raise IdentifierServiceError("identifiers.org resolver response has no payload")
    return resolution


__all__ = [
    "BARE_IDENTIFIER_NAMESPACE_ALLOWLIST",
    "DEFAULT_RESOLVER_TIMEOUT_SECONDS",
    "IDENTIFIERS_ORG_BASE_URL",
    "MAX_PROMPT_CANDIDATES",
    "RESOLVER_BASE_URL",
    "IdentifierNotFoundError",
    "IdentifierRef",
    "IdentifierResolutionError",
    "IdentifierServiceError",
    "IdentifierSyntaxError",
    "Resolution",
    "build_identifiers_org_url",
    "extract_identifier_candidates",
    "parse_identifier",
    "resolve_identifier",
]
