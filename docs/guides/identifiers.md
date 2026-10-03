# Identifiers.org over MCP

BioV exposes identifier-backed file descriptions, native provider metadata,
and the complete
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

Examples of data resource URIs include:

```text
refseq.gcf://GCF_000006945.2
uniprot://P12345
pubmed://22140103
clinvar://65533
dbsnp://rs121909098
geo://GSE1000
```

The four NCBI schemes validate accessions with their packaged identifiers.org
namespace rules. The GEO rule accepts GPL, GSM, GSE, and GDS accessions.

The RefSeq resource runs `datasets summary genome accession` and returns its
original genome-summary JSON stdout. It does not download a data package or
write the artifact cache; complete-package path resolution remains exclusive to
analysis code calling `biov.path`. The UniProt resource requests the complete
original UniProtKB REST JSON on a cache miss, stores those bytes under
`$BIOV_HOME/artifacts/uniprot/<accession>/`, and returns the cached content
unchanged. It does not request FASTA.

All other supported file namespaces, including PubMed, ClinVar, dbSNP, GEO,
PDB, Ensembl, ChEMBL, PubChem, arXiv, EMDB, ClinicalTrials.gov, DailyMed,
Reactome, and ENCODE, return a JSON description with `uri`, `default_kind`,
and `representations`. Reading that description validates the accession but
does not download the file or establish its availability. Use the URI with
`biov.path`, `biov.open`, or fsspec inside the analysis environment to obtain
the actual file. None of these resources returns an executor-local path.

For example, `pubmed://23193287` describes an `article_pdf`, while opening
that URI fetches the available PMC PDF. `geo://GSE100` describes the default
`expression_matrix` and explicit `soft` alternatives. Their actual download
limitations are documented in the [file provider guide](artifacts.md).

`clinvar://` accepts numeric Variation IDs; RCV and SCV accessions belong to
the separate `clinvar.record` and `clinvar.submission` registry namespaces and
continue to use `identifiers://` resolver resources.

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
packaged `identifiers_org_registry.json`; `biov update` is
the explicit command that refreshes that cached snapshot.

The RefSeq GCF registry rule is `^GCF_[0-9]{9}(\.[0-9]+)?$`. Consequently,
`GCF_000006945.2` is valid, while `GCF_00006945` is not because it has only
eight digits after `GCF_`. Network failures, malformed upstream responses, and
invalid accessions remain distinguishable errors. Other per-registry schemes
such as `go://...` and `doi://...`, and the former
`identifiers://resolve/...` resource do not exist.

## Parse IDs from a prompt

The `parse_identifiers` tool accepts the complete prompt as `prompt`. It finds
explicit Compact Identifiers, identifiers.org URLs, BioV resource URIs, and
explicitly allowlisted unambiguous bare IDs. It canonicalizes common variants,
validates each candidate against identifiers.org, preserves first-occurrence
order, removes duplicates, and returns links to the data scheme when supported
or otherwise to the generic identifiers resource. The first content block is
a JSON summary with `count` and `resource_uris`; the remaining blocks are
MCP resource links. A prompt with no valid identifiers returns a summary
with `count: 0` and an empty `resource_uris` array.

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
MCP embedded-resource content block. It accepts every supported identifier resource URI, including the file
namespaces listed by `biov.artifact_capabilities()`.

## Refresh the registry asset

The packaged `identifiers_org_registry.json` is the complete official resolver
response in its native nested shape. Runtime indexes read the required
namespace fields in memory without modifying the stored response.

Synchronize it with:

```console
biov update
```

The command validates the JSON envelope and fields needed by runtime routing,
then writes the upstream response bytes unchanged. A byte-identical response is
left untouched; a changed response atomically replaces the asset. Unknown
fields and nesting pass through unchanged. Use `--force` to rewrite an
unchanged response or `--output PATH` to target another asset file. Without
`--output`, the command writes the packaged asset of the running installation,
which is the checkout's `src/biov/assets/identifiers_org_registry.json` when
BioV is installed from this repository.

MCP resolution is discovery, not code execution. For analysis, keep the
identifier in generated source and resolve its path inside the target executor;
see [Identifier-backed analysis](artifacts.md).
