# Identifiers.org over MCP

BioV exposes two data-backed namespaces plus the complete
[identifiers.org](https://identifiers.org/) registry through MCP. The
integration uses official upstream data and does not require an API key.

## Start the server

After installing or syncing BioV, configure an MCP host to launch `biov mcp`
over standard input/output. For a checkout of this repository, a configuration
looks like this (replace the path):

```json
{
  "mcpServers": {
    "biov": {
      "command": "/absolute/path/to/biov/.venv/bin/biov",
      "args": ["mcp"]
    }
  }
}
```

The process reserves standard output for MCP messages.

## Read resources

The two data schemes return canonical provider data:

```text
refseq.gcf://GCF_000006945.2
uniprot://P12345
```

The RefSeq resource runs `datasets summary genome accession` and returns its
original genome-summary JSON stdout. It does not download a data package or
write the artifact cache; complete-package path resolution remains exclusive to
analysis code calling `biov.path`. The UniProt resource requests the complete
original UniProtKB REST JSON on a cache miss, stores those bytes under
`$BIOV_HOME/artifacts/uniprot/<accession>/`, and returns the cached content
unchanged. It does not request FASTA. Neither response contains a BioV wrapper
or an executor-local path.

All registry namespaces share the `identifiers` scheme:

```text
identifiers://go
identifiers://go:GO:0006915
identifiers://doi:10.1038/s41586-020-2649-2
identifiers://3dmet:B00162
```

`identifiers://<registry>` returns the complete native namespace object from
the packaged official registry snapshot, including provider resources and
institutions. `identifiers://<registry>:<id>` validates the local ID with that
namespace's rule and returns the complete native resolver JSON. It does not
fetch the resolved provider page or pretend that resolver metadata is the
entity's biological data.

Registry reads never synchronize over the network. They read the repository's
packaged `identifiers_org_registry.json`; `biov update-identifiers-registry` is
the explicit command that refreshes that cached snapshot.

The RefSeq GCF registry rule is `^GCF_[0-9]{9}(\.[0-9]+)?$`. Consequently,
`GCF_000006945.2` is valid, while `GCF_00006945` is not because it has only
eight digits after `GCF_`. Network failures, malformed upstream responses, and
invalid accessions remain distinguishable errors. Per-registry schemes such as
`go://...`, `doi://...`, and the former `identifiers://resolve/...` resource do
not exist.

## Parse IDs from a prompt

The `parse_identifiers` tool accepts the complete prompt as `prompt`. It finds
explicit Compact Identifiers, identifiers.org URLs, BioV resource URIs, and
explicitly allowlisted unambiguous bare IDs. It canonicalizes common variants,
validates each candidate against identifiers.org, preserves first-occurrence
order, removes duplicates, and returns links to the data scheme when supported
or otherwise to the generic identifiers resource.

Recognized forms include:

```text
uniprot:P12345
GO:0006915
ols/taxonomy:9606
doi:10.1038/s41586-020-2649-2
https://identifiers.org/pubmed:22140103
refseq.gcf://GCF_000006945.2
identifiers://doi:10.1038/s41586-020-2649-2
refseq.gcf:GCF_000006945.2
GCF_000006945.2
```

Compact Identifiers follow `[provider/]namespace:accession`. Provider-qualified
input is accepted but maps to a provider-independent resource URI. Explicit
resource URIs map back to the same Compact Identifier before validation.

Bare-ID inference is a curated namespace-prefix allowlist, not a scan across all
registry patterns. The current allowlist contains only `refseq.gcf`, whose
`GCF_` envelope is unambiguous. Thus bare `GCF_000006945.2` is recognized, but
bare `GCA_000155495.1` and `P12345` are not guessed. Equivalent URI, Compact,
and bare forms are resolved once. Syntactic candidates rejected by the official
resolver are omitted. Prompt work is capped at 50 unique candidates, and each
resolver request has a 10-second timeout.

## Tool-only clients

Some MCP clients can list and call tools but cannot call `resources/read`.
Pass a URI returned by `parse_identifiers` to:

```text
resolve_identifiers(uri="identifiers://go:GO:0006915")
```

The tool returns the same MIME type and content as the resource in a standard
MCP embedded-resource content block. It also accepts `refseq.gcf://` and
`uniprot://` resource URIs.

## Refresh the registry asset

The packaged `identifiers_org_registry.json` is the complete official resolver
response in its native nested shape. Runtime indexes read the required
namespace fields in memory without modifying the stored response.

Synchronize it with:

```console
biov update-identifiers-registry
```

The command validates the JSON envelope and fields needed by runtime routing,
then writes the upstream response bytes unchanged. A byte-identical response is
left untouched; a changed response atomically replaces the asset. Unknown
fields and nesting pass through unchanged. Use `--force` to rewrite an
unchanged response or `--output PATH` to target another asset file.

The version-controlled `scripts/update_identifiers_registry.py` delegates to
the same CLI command, defaulting to the checkout's asset.

MCP resolution is discovery, not code execution. For analysis, keep the
identifier in generated source and resolve its path inside the target executor;
see [Identifier-backed analysis](artifacts.md).
