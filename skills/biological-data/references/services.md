# Service-specific queries

These are task instructions and upstream documentation links. No local API
schema snapshot or BioV `query_database` wrapper is required.

## Sequence, structure and regulatory data

| Service | Query decision and native interface |
| --- | --- |
| UniProt | Search organism, reviewed status and accession explicitly with the [REST API](https://www.uniprot.org/help/api_queries); stream FASTA/TSV exports for batches. Preserve isoform IDs and requested columns. |
| AlphaFold | [Prediction metadata](https://alphafold.ebi.ac.uk/api-docs) identifies model files and confidence data. Require the intended model, sequence range and version when several models exist; predicted confidence is not experimental validation. |
| InterPro | Select protein-versus-entry queries using the [InterPro API](https://www.ebi.ac.uk/interpro/api/). A family ID identifies annotations and member sets, not one sequence. Follow pagination for member proteins and retain match coordinates. |
| PDB | Use the official [Search API](https://search.rcsb.org/) to find entry/entity IDs, then retrieve the chosen coordinate file. Distinguish entry, polymer entity, chain and assembly identifiers. |
| BLAST | Submit the provided sequence to [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html) or an explicitly selected local database with BLAST+. Record program, database version, E-value, scoring and hit limits; preserve complete native output and the job's actual completion state. |
| ENCODE | Search [ENCODE metadata](https://www.encodeproject.org/help/rest-api/) by assay, biosample, assembly and output type. Choose explicit ENCFF files before download; use the provider's href and MD5. For cCRE region queries specify assembly, interval convention and overlap criterion. |
| UCSC | Use the [Genome Browser API](https://genome.ucsc.edu/goldenPath/help/api.html) or published bigBed/bigWig files with assembly and track name. Record 0-based half-open coordinates and the exact track release. |
| Ensembl | Use the [REST documentation](https://rest.ensembl.org/) for lookup/overlap queries with species and assembly. Homology and multiple transcript results require explicit selection; do not silently choose the first transcript. |
| JASPAR | Select collection, species/taxon, profile version and PFM/PWM representation using the [official API](https://jaspar.elixir.no/api/). Matrix IDs are distinct from TF names; scanning also needs background frequencies and a stated threshold. |
| ReMap | Choose organism, genome build, release and experiment/peak aggregation from [official downloads](https://remap.univ-amu.fr/). Read the selected BED/peak file directly with fsspec; the service name alone does not identify a dataset. |
| RegulomeDB | Use the [official service](https://regulomedb.org/) with assembly-qualified variants/regions. Retain evidence categories and version; the score is regulatory evidence, not a causal conclusion. |
| EMDB | The [archive](https://www.ebi.ac.uk/emdb/documentation/faq) provides primary and half maps, metadata and validation reports. Use the named map type; never replace full-resolution data with a visualization mesh. |
| PRIDE | Use [pridepy](https://github.com/PRIDE-Archive/pridepy) or the [Archive API](https://www.ebi.ac.uk/pride/ws/archive/v2/docs/api-guide.html) to list all files for a PXD project and choose mzML/mzIdentML/mzTab or raw files deliberately. A project ID does not identify one file. |

## Variants, studies and clinical evidence

| Service | Query decision and native interface |
| --- | --- |
| ClinVar | Search with [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/) and retrieve full variation records. Preserve review status, condition, submitters, conflicts and dates. |
| dbSNP | Use [NCBI Variation Services](https://api.ncbi.nlm.nih.gov/variation/v0/) for RefSNP records and assembly-qualified placements. Merged/withdrawn records require explicit handling; an rsID can have multiple alleles. |
| GEO | [GEO downloads](https://www.ncbi.nlm.nih.gov/geo/info/download.html) separate Series matrices, SOFT and supplementary files. Match sample IDs and platform annotation; retain normalization metadata and do not label a raw-count file normalized. |
| GWAS Catalog | Use the [official API](https://www.ebi.ac.uk/gwas/rest/docs/api) for studies/associations and the summary-statistics file list for full results. Distinguish lead associations from genome-wide summary statistics; record ancestry, sample size and assembly. |
| gnomAD | Query the [native GraphQL API](https://gnomad.broadinstitute.org/api) with explicit dataset release and variant/gene ID. Preserve allele count/number, coverage, filters and population. Do not combine GRCh37 and GRCh38 releases without a defined conversion. |
| Open Targets | Use the [Platform API documentation](https://platform-docs.opentargets.org/data-access/graphql-api) and native GraphQL schema. Select target, disease and evidence source; association scores are not effect sizes or established causality. |
| cBioPortal | Select portal, study, sample list and molecular profile using its [OpenAPI documentation](https://www.cbioportal.org/api/swagger-ui/index.html). Join sample and patient IDs explicitly and retain missing values. |
| ClinicalTrials.gov | Use [API v2](https://clinicaltrials.gov/data-api/api) and its OpenAPI schema. Specify recruitment, condition and study type filters; follow pageToken. Keep protocol fields separate from reported results. |
| DailyMed | Use [SPL search and retrieval](https://dailymed.nlm.nih.gov/dailymed/webservices-help/v2/spls_api.cfm) with set ID and label version. Brand names can match several formulations. |
| openFDA | Use the [official API](https://open.fda.gov/apis/) for a named endpoint and query. Preserve endpoint-specific identifiers and reporting limitations; spontaneous adverse-event counts are not incidence rates. |
| Synapse | Use [synapseclient](https://python-docs.synapse.org/) with an explicit entity/version and the user's configured credentials. Respect dataset access requirements; keep credentials out of results and never substitute a public dataset for denied data. |
| MPD | Choose phenotype, strain, sex and measurement units from the [Mouse Phenome Database](https://phenome.jax.org/). Preserve study identifiers and protocol differences when combining results. |

## Molecules, pathways and organism resources

| Service | Query decision and native interface |
| --- | --- |
| ChEMBL | Use the [official client/API](https://chembl.gitbook.io/chembl-interface-documentation/web-services) for molecules, targets, assays or activities. Keep assay type, target confidence, relation operators and units; an IC50 is not interchangeable with Ki. |
| PubChem | Use [PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) with CID/SID/AID chosen deliberately. Structure download is separate from assay results and source attribution; retain stereochemistry and salt/parent identity. |
| UniChem | Use [UniChem](https://www.ebi.ac.uk/unichem/) for cross-source structure identifiers. A mapping is a source-specific identifier relationship, not a potency or biological-equivalence result. |
| Guide to Pharmacology | Use the [official web services](https://www.guidetopharmacology.org/webServices.jsp) for targets/ligands/interactions, recording species and affinity type/units. |
| KEGG | Use the [documented API](https://www.kegg.jp/kegg/rest/) for explicit organism-qualified genes, pathways, compounds and their FASTA/KGML/mol outputs. The public API is restricted to academic use; non-academic use requires the provider's agreement. Do not infer permission from Biomni's old configuration. |
| Reactome | Use the [Content Service](https://reactome.org/dev/content-service) for events, species and SBML exports. Enrichment is a separate analysis requiring the submitted universe and method; preserve the actual service response and release. |
| QuickGO | Use the [official API](https://www.ebi.ac.uk/QuickGO/api/index.html) for ontology terms or annotations. Specify taxon, evidence codes, gene-product namespace and release; follow pagination/export limits. |
| STRING | Use the [official API](https://string-db.org/help/api/) with species and version. Record whether edges are physical or functional, score channels and threshold; enrichment requires the tested background. |
| Monarch | Use the [official API](https://api.monarchinitiative.org/docs) with typed phenotype/disease/gene identifiers. Preserve relation and evidence provenance; a graph path alone is not a causal mechanism. |
| IUCN | Use the [current API](https://api.iucnredlist.org/) with authorized credentials, taxon and assessment scope/year. National/regional and global assessments are distinct. |
| Paleobiology | Use the [PBDB Data Service](https://paleobiodb.org/data1.2/) with taxon, time interval, geography and requested fields. Preserve collection context and age uncertainty. |
| WoRMS | Use the [REST service](https://www.marinespecies.org/rest/) to resolve names to AphiaIDs and accepted taxa. Preserve synonym/accepted-name status, authority and marine/freshwater scope. |

## Literature

Search PubMed with [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
and arXiv with its [native API](https://info.arxiv.org/help/api/user-manual.html),
recording the query, date, pagination and exact IDs. General web/Scholar searches
belong to the calling agent's search tools. Retrieve supplementary data through
the article/publisher's named links and retain the parent DOI, supplement label,
file format and license. A PDF, abstract and supplementary dataset are distinct
artifacts; absence of one must not be concealed by returning another.
