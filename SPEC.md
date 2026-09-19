# SPEC

## §G GOAL
Own interval/sequence APIs; expose persistent-ID discovery, environment-local biological artifacts & whole-script execution to LLMs without wrapping analysis libraries.

## §C CONSTRAINTS
- Python `>=3.12`; RuRanges only interval kernel; Biopython only sequence-algorithm backend
- interval semantics ≠ sequence semantics; ⊥ new DataFrame subclass or large framework
- ⊥ PyRanges objects, conversion paths, compatibility aliases, or legacy interval dependencies
- public behavior ! documented & acceptance-tested
- identifiers.org input = explicit Compact Identifier, identifiers.org URI, or BioV resource URI; bare ID inference ! curated unambiguous namespace allowlist; ⊥ registry-wide schema stripping
- identifiers.org resolution access ! fixed to `https://resolver.api.identifiers.org`; MCP transport = stdio
- MCP schemes = data-backed `refseq.gcf|uniprot` + generic `identifiers`; ⊥ per-registry schemes or `identifiers://resolve/{+compact_id}` compatibility
- BioV integration scales by identifier namespace × artifact kind; ⊥ wrappers for Biopython/other analysis functions
- generated analysis ! use ordinary ecosystem APIs; complete script executes in selected local|LSF environment so concrete storage paths never cross executor boundary
- RefSeq assembly path resolution ! official `datasets download genome accession` CLI with `--include gff3,rna,cds,protein,genome,seq-report`; extracted data package structure ! preserved verbatim
- RefSeq MCP metadata ! official `datasets summary genome accession` JSON stdout; ⊥ artifact path resolution, package download/extraction or cache write
- UniProt path resolution ! official full-entry `.json` and protein `.fasta` responses cached independently on demand; bytes & accession filenames ! preserved; ⊥ field-filtered metadata or implicit PDB/AlphaFold selection

## §I INTERFACES
- api: `BioDataFrame.overlap(other, how, seqid_col, start_col, end_col, strand_col)` → selected self rows
- api: `BioDataFrame.intersect(other, seqid_col, start_col, end_col, strand_col)` → clipped self rows per overlap pair
- api: `BioDataFrame.subtract_ranges(other, seqid_col, start_col, end_col, strand_col)` → residual self fragments
- api: `BioDataFrame.nearest(other, seqid_col, start_col, end_col, strand_col, suffix, how)` → self + nearest other columns + `Distance`
- dtype: `biov.dna` | `biov.rna` | `biov.protein` → validated nullable uppercase sequence storage
- accessor: `Series.seq.length` → nullable integer Series
- accessor: DNA/RNA `Series.seq.reverse_complement()` | `gc_fraction()` | `translate(table=1, to_stop=False)`
- accessor: protein `Series.seq.molecular_weight()` | `isoelectric_point()` → nullable float Series; `amino_acid_composition()` → 20-column percentage DataFrame
- resource: `refseq.gcf://{+accession}` → original NCBI genome-summary JSON; `uniprot://{+accession}` → original complete UniProtKB JSON
- resource: `identifiers://{registry}` → native registry namespace object; `identifiers://{registry}:{+id}` → native resolver JSON
- tool: `parse_identifiers(prompt)` → ordered unique MCP resource links for resolver-valid explicit IDs, BioV resource URIs & allowlisted bare IDs
- tool: `resolve_identifiers(uri)` → the same resource content as a standard MCP embedded resource for tool-only clients
- cmd: `biov mcp` → BioV MCP server over stdio; ⊥ standalone `biov-mcp`
- cmd: `biov update-identifiers-registry` → fetch & atomically refresh the packaged registry asset
- file: `src/biov/assets/identifiers_org_registry.json` → original complete resolver-dataset JSON response body
- script: `scripts/update_identifiers_registry.py` → fetch, validate & atomically replace raw registry asset
- api: `parse_identifier(value)` → exact single `IdentifierRef`; same generated URI, Compact ID, identifiers.org URL & curated bare-ID syntax as prompt parser
- api: `path(identifier, artifact=None)` → namespace-default environment-local immutable `Artifact` implementing `os.PathLike[str]`
- api: `open(identifier, artifact=None, mode="rb")` → handle for cached artifact
- provider: `refseq.gcf` × `genome_fasta` → original catalog-selected genomic FASTA inside a complete NCBI Datasets package
- provider: `refseq.gcf` × `annotation_gff3` → original catalog-selected GFF3 inside a complete NCBI Datasets package
- provider: `refseq.gcf` × `rna_fasta` → original catalog-selected RNA FASTA inside a complete NCBI Datasets package
- provider: `refseq.gcf` × `cds_fasta` → original catalog-selected CDS FASTA inside a complete NCBI Datasets package
- provider: `refseq.gcf` × `protein_fasta` → original catalog-selected protein FASTA inside a complete NCBI Datasets package
- provider: `uniprot` × `protein_fasta` → original accession-named UniProtKB FASTA response
- provider: `uniprot` × `entry_json` → original complete accession-named UniProtKB JSON response
- cmd: `biov run [OPTIONS] SCRIPT [ARGS]...` → run complete Python script locally or submit it to LSF
- env: `BIOV_LSF_PYTHON` ? Python executable visible from LSF execution hosts; default = submitting interpreter
- file: `$BIOV_HOME/artifacts/refseq.gcf/<requested_accession>/` → unmodified extracted NCBI Datasets package root (`README.md`, `md5sum.txt`, `ncbi_dataset/...`)
- file: `$BIOV_HOME/artifacts/uniprot/<accession>/` → independently cached unmodified `<accession>.fasta` and/or complete `<accession>.json`; ⊥ manifest
- file: `.github/workflows/ci.yml` → test matrix, hooks & wheel/sdist inspection on push and pull requests
- file: `.github/workflows/registry-drift.yml` → scheduled asset sync; opens a pull request when upstream changes

## §R RESEARCH
id|topic|finding|src
R1|RuRanges surface|stateless `ruranges.numpy` functions accept/return NumPy arrays; groups integer-coded; strand boolean only where required|https://github.com/pyranges/ruranges_py
R2|RuRanges kernels|`overlaps`, `nearest`, `subtract` return source indices; nearest physical directions = `forward`/`backward` & overlap distance = 0|https://raw.githubusercontent.com/pyranges/ruranges_py/master/ruranges/numpy.py
R3|pandas extension|custom dtype + 1-D ExtensionArray preserve semantic type; accessor init ! reject wrong dtype with `AttributeError`|https://pandas.pydata.org/docs/development/extending.html
R4|Biopython sequence math|weighted GC defines ambiguous IUPAC handling & empty GC = 0; molecular weight requires unambiguous residues|https://biopython.org/docs/latest/api/Bio.SeqUtils.html
R5|identifiers.org resolution|Compact Identifier = `[provider/]namespace:accession`; resolver `GET https://resolver.api.identifiers.org/{COMPACT_ID}` returns provider resources|https://docs.identifiers.org/pages/api.html
R6|MCP Python SDK|v2 `MCPServer` exposes typed tools & URI-template resources over stdio|https://py.sdk.modelcontextprotocol.io/
R7|identifiers.org registry dataset|`GET /resolutionApi/getResolverDataset` returns complete registry including namespace patterns, provider resources & institutions|https://docs.identifiers.org/pages/api.html
R8|NCBI genome CLI|`datasets download genome accession <GCF> --include ... --filename <zip> --no-progressbar` downloads one official ZIP; include values cover genome, GFF3, RNA, CDS, protein & sequence report|https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/
R9|NCBI genome package|extracted package root contains `README.md`, `md5sum.txt`, `ncbi_dataset/data/...`; each assembly keeps original files under its accession directory & catalog labels genome FASTA `GENOMIC_NUCLEOTIDE_FASTA`|https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/
R10|LSF submission|`bsub` acceptance assigns job ID while job may remain `PEND`; `DONE|EXIT` occur later ∴ ordinary submission ≠ completion|https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=management-job-lifecycle
R11|UniProt individual entry|`GET https://rest.uniprot.org/uniprotkb/<accession>.fasta` is the documented direct retrieval form for one UniProtKB entry in FASTA|https://www.uniprot.org/help/api_retrieve_entries
R12|UniProt structure links|one UniProtKB entry can expose many PDB cross-references with distinct methods, resolutions & chain coverage ∴ UniProt accession ≠ unique PDB coordinate file|https://rest.uniprot.org/uniprotkb/P42212.json?fields=xref_pdb
R13|UniProt complete entry|individual-entry REST retrieval supports accession-qualified `.json`; unfiltered response retains full entry metadata including database cross-references|https://www.uniprot.org/help/api_retrieve_entries
R14|NCBI genome summary|`datasets summary genome accession <GCF>` returns assembled-genome metadata as JSON without downloading a genome data package|https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/summary/genome/datasets_summary_genome_accession/

## §V INVARIANTS
V1: intervals use 0-based, end-exclusive `[start,end)`; integer `0 ≤ start < end`; touching boundaries ≠ overlap
V2: operations group by exact `seqid`; when `strand_col` exists on both, only `+|-` valid & group by exact `seqid,strand`; one-sided/invalid strand → `ValueError`; `strand_col=None` ignores strand
V3: outputs preserve self input order; pair expansions use self row then other row order; fragments use self row then ascending coordinate; duplicate input rows remain distinct
V4: empty self/other inputs return typed empty or unchanged results with stable schema; ⊥ kernel panic
V5: `overlap` returns each matching self row once; `first|last` select membership, `containment` means self contains other, `member` means self contained by other
V6: `intersect` emits one clipped self-metadata row per overlap pair; overlap duplicates yield duplicate output rows
V7: `subtract_ranges` emits every non-empty residual fragment with self metadata; overlapping/duplicate masks do not duplicate residual space
V8: `nearest` returns ≤1 row/query; overlap `Distance=0`, otherwise `gap+1` so adjacent half-open intervals have `Distance=1`; `next|previous` = genomic right|left; `upstream|downstream` = strand-aware 5′|3′; tie → lowest other input row; missing group candidate → omit query
V9: public package & lock contain ⊥ `pyranges`, `sorted-nearest`, `ncls`; unnecessary `setuptools` absent
V10: ordinary string Series `.seq` → `AttributeError`; caller ! choose/carry `biov.dna|rna|protein`; ⊥ content guessing
V11: sequence storage uppercases valid strings, preserves `pd.NA`, accepts empty strings & declared IUPAC alphabets, rejects non-string/invalid symbols with `SequenceValidationError`
V12: DNA/RNA reverse complement preserves dtype/nulls; weighted GC handles IUPAC & empty string; translation requires complete codons, honors `table,to_stop`, preserves nulls, returns `biov.protein`
V13: protein mass/pI/composition use non-empty canonical 20 amino acids; extended IUPAC, stop, or empty sequence stored but analysis → stable `SequenceValidationError`; null results remain null; composition columns = canonical amino-acid order & values = percentages
V14: RuRanges called only from interval module; Biopython algorithms called only from sequence module
V15: pre-existing public BioV behavior & tests remain intact
V16: prompt parser recognizes explicit Compact Identifiers, identifiers.org URLs, BioV resource URIs & allowlisted bare IDs; canonicalizes variants, preserves first occurrence order, deduplicates, accepts provider/slash accessions & trims prose punctuation
V17: parser emits links only for resolver-valid IDs; invalid candidates omitted; upstream/service/JSON failures surface distinctly ≠ invalid ID
V18: `identifiers://registry:id` round-trips reserved accession characters & returns native resolver JSON; `identifiers://registry` losslessly reconstructs the native namespace object; invalid registry/accession → stable resource error
V19: identifiers.org integration is read-only/idempotent, accesses only fixed resolver origin, has bounded prompt candidates & request timeout
V20: `biov mcp` starts stdio without protocol-corrupting stdout; `biov-mcp` ∉ installed scripts; existing public APIs remain importable
V21: asset preserves complete official response in native nested shape; namespace resources, institutions & locations remain upstream-owned objects
V22: MCP publishes exactly four templates across three schemes: `refseq.gcf`, `uniprot`, and the two `identifiers` forms; per-registry and old `identifiers://resolve/*` templates ∉ resources
V23: data template metadata exposes prefix, accession regex, sample & namespace ID; every ID read validates the asset regex; `namespaceEmbeddedInLui` reconstructs resolver input correctly
V24: prompt tool uses resolver-parsed namespace/local ID → data resource for `refseq.gcf|uniprot`, otherwise generic identifiers resource; provider-qualified input maps to a location-independent URI; output URI deduplicated in prompt order
V25: asset validation checks only fields required for runtime indexing/routing; unknown upstream fields & nesting ! preserved; asset ! packaged in wheel/sdist
V26: registry update validates upstream JSON/runtime fields; byte-identical body skips replacement; changed body atomically replaces asset unchanged
V27: bare-ID recognition checks only an explicit namespace-prefix allowlist; `refseq.gcf` ∈ allowlist; non-allowlisted registry patterns ∉ inference; invalid shorthand omitted before resolver access
V28: each accepted identifier syntax variant ! pass an isolated acceptance case without another valid fallback form; combined variants ! canonicalize & deduplicate before resolver access
V29: `parse_identifier` accepts exactly one full reference, validates packaged namespace regex, returns direct data URI for supported namespaces or generic identifiers URI otherwise without resolver I/O; prose, multiple references, per-registry schemes, unknown namespaces & non-allowlisted bare IDs → stable syntax error
V30: artifact support registry keyed only by `(namespace,artifact kind)`; analysis remains ordinary Biopython/ecosystem code; ⊥ BioV GC/alignment/etc wrapper proliferation
V31: versioned GCF path request returns that exact catalog assembly; versionless GCF accepts exactly one matching versioned catalog assembly; missing/ambiguous/non-GCF package → stable artifact error
V32: `Artifact` is accepted wherever `os.PathLike[str]` is accepted & exposes original package member path, package root, requested/canonical identifier, kind & byte size
V33: valid package-directory cache hit reads only local catalog/file metadata & skips `datasets`; cache miss downloads/extracts in a unique sibling staging directory then atomically publishes the complete package root; failed/interrupted writes never become cache hits
V34: command argv = `datasets download genome accession <accession> --include gff3,rna,cds,protein,genome,seq-report --filename <temporary-zip> --no-progressbar`; every safe ZIP member retains its exact relative path; ⊥ renamed/copied FASTA, BioV manifest, partial extraction or direct REST download
V35: original package catalog selects exactly one existing `GENOMIC_NUCLEOTIDE_FASTA`; unsafe/duplicate ZIP paths, invalid catalog, accession mismatch or missing/duplicate FASTA → stable package error; missing/rejected `datasets` CLI → stable service error
V36: `biov run` passes Python executable, absolute script & arguments as argv without shell interpolation; complete script—including `path`—runs inside selected executor
V37: local executor inherits current cwd/environment/stdio & command exit code becomes `biov run` exit code
V38: LSF executor submits via `bsub`, pins cwd, supports queue/name/stdout/stderr + executor-visible Python, returns parsed numeric job-ID receipt; accepted submission never claims job completion; missing/rejected/unknown receipt → stable execution error
V39: cached path is valid only inside current executor; LSF submission does not resolve or return compute-node paths to submitter; shared env/cache/network availability remains deployment configuration
V40: unsupported namespace × artifact, invalid identifier, assembly mismatch, malformed package, missing/rejected `datasets`, missing executor & rejected submission have distinct public exception types/messages
V41: omitted artifact kind dispatches by namespace: `refseq.gcf` → `genome_fasta`; `uniprot` → `protein_fasta`; explicit unsupported pair remains a stable `UnsupportedArtifactError`
V42: UniProt request URLs = fixed HTTPS origin + `/uniprotkb/<percent-encoded-accession>.fasta|.json`; streamed response bytes remain unchanged at accession-named files; ⊥ generated manifest, sequence/JSON rewrite or field filtering
V43: each UniProt artifact validates only its local regular file, bounded `sp|tr` FASTA header or JSON object, and exact accession; a miss stages and atomically publishes only the requested representation; absent/invalid/mismatched response never becomes a cache hit
V44: `uniprot://<accession>` defaults to sequence path while its cache retains complete entry metadata; multiple PDB cross-references or one predicted model ! require separately identified coordinate artifacts; ⊥ arbitrary structure selection
V45: `entry_json` and `protein_fasta` cache independently; requesting one never fetches or requires the other; full raw JSON retains all upstream PDB IDs and cross-reference properties
V46: registry asset bytes = upstream response body; ⊥ wrapper, flattening, foreign keys or derived fields; runtime indexes native records in memory without mutating them
V47: `refseq.gcf://<accession>` executes only `datasets summary genome accession <accession>`; validates one matching report then returns stdout unchanged; ⊥ `path`, package download/extraction or `$BIOV_HOME` write; `path(refseq.gcf)` remains V31–V35
V48: each registered `refseq.gcf` artifact kind maps to exactly one official catalog `fileType`; the catalog selects exactly one existing member per request; missing/duplicate members → stable package error; one cached package serves every kind

## §T TASKS
id|status|task|cites
T1|x|write contracts & failing acceptance tests|I.*,V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13
T2|x|replace interval adapter with RuRanges NumPy kernels|I.overlap,I.intersect,I.subtract_ranges,I.nearest,V1,V2,V3,V4,V5,V6,V7,V8,V9,V14
T3|x|add typed sequence EA/dtypes/accessor|I.dtype,I.accessor,V10,V11,V12,V13,V14
T4|x|update public docs, exports, dependencies & lock|V9,V15
T5|x|run full verification matrix|V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15
T6|x|write identifiers.org MCP acceptance tests|I.resource,I.tool,I.cmd,V16,V17,V18,V19,V20
T7|x|implement resolver client, resource, parser tool & stdio server|I.resource,I.tool,I.cmd,V16,V17,V18,V19,V20
T8|x|document MCP setup, parsing scope & error behavior|I.resource,I.tool,I.cmd,V16,V17,V18,V19,V20
T9|x|run full verification matrix|V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20
T10|x|write raw registry asset & generated-resource acceptance tests|I.resource,I.file,I.script,V18,V21,V22,V23,V24,V25
T11|x|persist raw registry response & generate full asset|I.file,I.script,V21,V22,V25
T12|x|generate namespace MCP resources & prompt links from asset; remove old URI|I.resource,I.tool,V16,V17,V18,V19,V22,V23,V24
T13|x|document registry URI schemes, aliases & asset refresh|I.resource,I.file,I.script,V22,V23,V24,V25
T14|x|run full verification matrix & package inspection|V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20,V21,V22,V23,V24,V25
T15|x|add `biov update-identifiers-registry` command & acceptance tests|I.cmd,I.file,V21,V25,V26
T16|x|replace standalone `biov-mcp` with `biov mcp` & update docs/package|I.cmd,V15,V20
T17|x|support generated URI variants & allowlisted unambiguous bare IDs in `parse_identifiers`|I.tool,V16,V17,V19,V24,V27,V28
T18|x|write single-ID, path resolution/cache/package & local/LSF execution acceptance tests|I.api,I.provider,I.cmd,V29,V30,V31,V32,V33,V34,V35,V36,V37,V38,V39,V40
T19|x|implement `IdentifierRef`, artifact provider registry, initial NCBI GCF FASTA provider/cache & public exports|I.api,I.provider,I.file,V29,V30,V31,V32,V33,V34,V35,V39,V40
T20|x|implement whole-script local/LSF executors & `biov run`|I.cmd,I.env,V36,V37,V38,V39,V40
T21|x|document LLM/Biopython usage, executor boundary, artifact cache & deployment requirements|I.api,I.cmd,I.env,V30,V31,V32,V38,V39,V40
T22|x|run full tests, hooks, package inspection & bounded official NCBI smoke validation|V15,V29,V30,V31,V32,V33,V34,V35,V36,V37,V38,V39,V40
T23|x|replace REST/flattened-FASTA acceptance tests with official CLI argv, complete-package layout, catalog selection & cache tests|I.api,I.provider,I.file,V31,V32,V33,V34,V35
T24|x|replace RefSeq REST adapter/manifest cache with `datasets` CLI & verbatim package cache|I.api,I.provider,I.file,V31,V32,V33,V34,V35
T25|x|document CLI prerequisite/package layout & run full verification|V15,V31,V32,V33,V34,V35
T26|x|write UniProt URI, official endpoint, raw-file cache & failure acceptance tests|I.api,I.provider,I.file,V27,V30,V41,V42,V43,V44
T27|x|implement namespace-default artifact dispatch & UniProt FASTA provider|I.api,I.provider,I.file,V30,V41,V42,V43,V44
T28|x|document UniProt/PDB boundary, run full verification & download official GFP `P42212`|V15,V41,V42,V43,V44
T29|x|backprop FASTA-only cache bug; preserve & expose complete UniProt JSON beside FASTA|I.api,I.provider,I.file,V42,V43,V44,V45
T30|x|replace relational registry snapshot with raw upstream response|I.file,I.script,V21,V25,V26,V46
T31|x|decouple RefSeq MCP summary reads from complete analysis-package path resolution|I.resource,V23,V47
T32|x|harden review findings: default MCP resource security, mapped registry CLI errors, bounded datasets/bsub subprocesses|V19,V20,V35,V38,V40
T33|x|generalize RefSeq catalog selection to annotation/rna/cds/protein artifact kinds|I.provider,V31,V32,V33,V35,V48
T34|x|add CI test/hook/package workflows & scheduled registry-drift sync; restore byte-exact asset|V15,V25,V26,V46

## §B BUGS
id|date|cause|fix
B1|2026-08-30|mixed-form acceptance prompt contained Compact ID fallback, masking absent URI/bare parsing|V28
B2|2026-08-30|new allowlist error path exposed incomplete MCP helper docstring contracts|docstrings completed; no behavior invariant
B3|2026-08-31|RefSeq provider reimplemented NCBI REST and flattened one renamed FASTA instead of retaining the official CLI data package|V31,V32,V33,V34,V35
B4|2026-09-01|UniProt provider conflated default computation artifact with complete upstream cache & downloaded FASTA only|V45
B5|2026-09-03|registry snapshot transformed upstream JSON into an unused relational format|V46
B6|2026-09-03|RefSeq MCP metadata read reused `genome_fasta` path resolution & downloaded complete analysis package|V47
