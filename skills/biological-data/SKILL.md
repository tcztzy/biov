---
name: biological-data
description: Find and retrieve biological sequences, structures, variants, expression data, papers, and study records. Use for selecting a database, resolving an accession, obtaining analysis files, or constructing an explicit official API query.
---

# Retrieve biological data

For software discovery and script execution, use
[scientific-software](../scientific-software/SKILL.md). Its catalog includes
the upstream software descriptions and links to BioV's execution instructions.

Distinguish discovery from retrieval. A gene symbol, disease, or natural-language
question may match many records. Search the relevant service with the organism,
assembly, dataset release, assay, or cohort needed by the question, then retain
the exact returned identifiers. A matching identifier regex does not establish
that a record exists or that it is the intended biological entity.

## Known identifiers

Use BioV in the execution environment to obtain files. `biov.path(uri)` returns
a path-like artifact whose `.path` is the actual local filename;
`biov.open(uri)` and `fsspec.open(uri)` return readable file objects.
`biov.artifact_capabilities()` lists the supported file kinds and defaults.
The existing `read_fasta` returns sequence records, not a DataFrame.

```python
import biov
import fsspec
from Bio.PDB import MMCIFParser

coordinates = biov.path("pdb://1CRN")
structure = MMCIFParser(QUIET=True).get_structure("1CRN", coordinates.path)
with fsspec.open("pubmed://23193287", "rb") as article:
    pdf_header = article.read(5)
```

| Identifier | File available through BioV/fsspec | Selection constraint |
| --- | --- | --- |
| `refseq.gcf://GCF_000006945.2` | Genomic FASTA; GFF3, RNA/CDS/protein FASTA via `artifact=` | Exact version when given; retain the complete original NCBI package. |
| `uniprot://P42212` | Protein FASTA; `entry_json`; `alphafold_cif` | AlphaFold is an explicit predicted-coordinate request and requires one matching model. Do not select one of many PDB cross-references arbitrarily. |
| `pubmed://23193287` | PDF from the current available PMC version | Some papers have no downloadable PDF. Report absence; do not substitute a summary. |
| `arxiv://1706.03762` | PDF | Include the version in the ID when reproducibility needs a fixed revision. |
| `clinvar://9` | Complete VCV XML | Evidence/assertions, not an automatic clinical interpretation. |
| `dbsnp://rs328` | Complete RefSNP JSON | Check assembly and allele placement before joining to genomic data. |
| `geo://GSE100` | Original Series Matrix table; explicit `soft` for GSE/GDS/GPL/GSM | Default requires one nonempty matrix. Multiple matrices need explicit source-file selection; no matrix is guessed. |
| `pdb://1CRN` | Original mmCIF | Asymmetric unit coordinates; do not call it a selected biological assembly. |
| `emdb://EMD-1001` | Decompressed primary density map | Preserve voxel size, origin and axis conventions. |
| `ensembl://ENSG00000139618` | Native FASTA | Ensembl's default sequence type depends on feature type; inspect the header. |
| `chembl.compound://CHEMBL25`, `pubchem.compound://2244` | Molecular structure file | PubChem uses 2D coordinates. Structure retrieval does not retrieve all bioactivity measurements. |
| `clinicaltrials://NCT00222573` | Complete study JSON | Distinguish registered protocol from posted results. |
| `dailymed://1efe378e-fee1-4ae9-8ea5-0fe2265fe2d8` | Complete SPL XML | Match the specific product/set ID. |
| `reactome://R-HSA-199420` | SBML export | Requires an exportable event; an SBML graph does not establish kinetic parameters. |
| `encode://ENCFF002CTW` | Original file with published name and compression | Requires an ENCFF file ID. An experiment ID does not identify one file. |

MCP `parse_identifiers` discovers identifiers and `resolve_identifiers` embeds
resource results for tool-only clients. New file resources describe their URI
and available representations as JSON; obtaining the actual file uses the
Python/fsspec interfaces above. RefSeq and UniProt MCP resources retain their
native provider metadata. Do not treat file-description JSON as the data file.

## Searches and services requiring more context

Use the [service-specific instructions](references/services.md) for searches,
cohort queries, genomic regions, batch operations, restricted datasets and
services without a unique file per accession. Use their maintained SDK or
official HTTP/GraphQL interface directly. Read the provider's current OpenAPI
description or GraphQL schema when supplied; BioV does not maintain a second
API catalog or a generic query dispatcher.

Follow the service's actual pagination or asynchronous-job protocol. A first
page is not a complete dataset. Request the fields needed for the analysis;
retain original output and request parameters before any scientific filtering.
Credential errors and unavailable records are failures, not successful result
strings. Do not automate real laboratory devices as part of data retrieval.

Return a JSON summary identifying the request, selected IDs, service/release,
file paths and formats, transformations, and actual missing/ambiguous results.
Keep full sequences, matrices, PDFs and coordinate files as files. The caller
chooses the scientific interpretation and downstream analysis.

Migration source: Biomni commit `400c1f366b96a35ca253e13c9b06c5076af41d65`,
`biomni/tool/database.py` and literature retrieval functions. This skill
replaces their embedded LLM query selection, not their historical signatures.
