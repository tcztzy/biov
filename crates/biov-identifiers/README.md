# biov-identifiers

A small offline Rust contract for biological references. This crate is independent
of Python, MCP, files, and network access. It does **not** replace Python BioV's
full identifiers.org registry/parser.

## Supported subset

- `refseq.gcf`: `GCF_` followed by exactly nine ASCII digits, optionally a dot and
  one or more ASCII digits, matching the packaged identifiers.org GCF rule.
  Revision text is preserved, including leading zeroes, without numeric overflow.
  An unversioned accession is not a pinned assembly revision. Passing local
  syntax validation does not show that an accession or revision exists.
- `uniprot`: the provider's six- and ten-character UniProtKB accession patterns.
  The parser cannot distinguish primary from secondary accessions or establish
  whether an entry is current. Isoforms (`-2`), dot revisions, processed-protein
  fragments (`#PRO_...`), and entry names are deliberately outside this subset.
  UniProt sequence/entry revisions require separate provider metadata.

`IdentifierRef::parse` / `FromStr` accepts compact IDs, BioV biological resource
URIs, generic `identifiers://namespace:accession` URIs, and identifiers.org
HTTP(S) URLs using the compact or legacy slash path. Only GCF is allowlisted
for bare accessions, matching Python BioV. Namespace, URL scheme, and URL host
comparison is case-insensitive; accession case is preserved and checked.

Only one complete reference is accepted. Outer whitespace is trimmed. Prose,
multiple IDs, provider-prefixed IDs, percent-escaped paths, query strings,
fragments, trailing slashes, and unsupported namespaces are rejected. The
supported accessions need no percent escaping. No suffix or URL qualifier is
silently removed. This is a deliberately narrower contract than the Python
parser, especially for the broad, imperfect upstream UniProt registry regex.

## API and data integrity

`Namespace::{RefSeqGcf, UniProt}` and `IdentifierRef::new(namespace, accession)`
allow callers to provide an explicit namespace. Read-only accessors expose the
namespace, full accession, base accession, and optional GCF revision. Compact
IDs and resource URIs are derived from those coordinates. Equality includes
namespace and the complete accession, so unversioned/versioned or differently
versioned references are not collapsed.

Serde stores only `namespace` and `accession`; deserialization validates them
and rejects unknown fields. There is no independently writable version or URI
that could disagree with the accession. The only runtime dependency is serde.

A biological reference identifies an external biological record. An opaque
BioV dataset handle identifies a managed local result; these are not
interchangeable. Dataset handles and paths are rejected. Consumers should store
the validated biological reference as source/provenance metadata alongside,
not in place of, their separately typed dataset handle.

## Scientific contract and sources

Tests adapt the existing Python cases in `tests/test_artifacts.py` and
`tests/test_identifiers_mcp.py`, the packaged GCF registry rule, and examples
from the primary provider documentation:

- [NCBI assembly versioning](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genome-processing/version-status/)
- [UniProt accession syntax](https://www.uniprot.org/help/accession_numbers)

The tests verify syntax and preservation of version meaning. They do not claim
record existence, resolution, annotation version, sequence identity, or current
provider status. In particular, GCF and GCA assemblies are distinct and may
have different revision numbers; this crate does not translate between them.

From the repository workspace, run `cargo test -p biov-identifiers`.
