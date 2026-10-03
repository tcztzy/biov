---
name: analyze-omics
description: Analyze bulk genomic, epigenomic, variant, association, and chromatin data with established tools. Use for peak calling, differential accessibility, enrichment, fine-mapping, copy number, somatic variants, orthology, or genomic prediction.
---

# Genomic and epigenomic analysis

Require sample metadata, reference build/version, biological replicates and
the tested contrast before choosing a method. Match sample and feature IDs
explicitly. Use native tool outputs (VCF, BED, cool/mcool, H5AD, TSV/Parquet)
plus JSON containing input references, parameters, software versions and QC.
Run the existing scientific tools directly, not the retired Biomni wrappers.

## Variant and association analysis

| Task | Inputs and method | Checks and outputs |
| --- | --- | --- |
| Coordinate liftover | Source and destination assembly, interval convention and matching UCSC chain; use CrossMap. | Mapped, unmapped and multiply mapped intervals with strand. Do not silently select the first mapping or retain the original coordinates on failure. |
| Somatic small variants | Matched tumor/normal BAM/CRAM and matching reference; use GATK Mutect2 with the relevant germline resources, contamination/orientation modelling and FilterMutectCalls. Annotate using the selected genome/transcript release. | Filtered annotated VCF and QC, including contamination and callable coverage. Tumor-only analysis needs its own stated assumptions; do not improvise a normal sample. |
| Structural variants | LUMPY is available for appropriate paired-end short-read data; other read technologies require a separately selected caller. Annotate the resulting VCF with a specified annotation tool and versioned sources. | SV types, breakpoints, support, filters and annotation source. Do not invoke the upstream placeholder `annotate_sv.py`. |
| Copy number | Use CNVkit with an appropriate reference, panel targets/antitargets or WGS mode; retain bins, segments and normalization QC. | Log2 copy-ratio segments and focal overlaps. Absolute copy number, purity, ploidy, LOH and HRD need models and data that identify them; log2 thresholds or segment counts alone do not. |
| Fine-mapping | Align GWAS effect alleles, variant order, ancestry, sample size and signed LD; use `susieR::susie_rss` or another justified fine-mapper. | PIPs, per-effect credible sets, coverage/purity and convergence diagnostics. Do not reuse the upstream Bernoulli-sampling neural objective or truncate a global PIP sum below the requested coverage. |
| Genomic prediction | Genotypes, phenotypes, covariates and a stated additive/dominance relationship model; use a maintained mixed-model implementation such as sommer. | Estimated variance components and predictions with family/population-aware held-out validation. In-sample phenotype correlation is not predictive accuracy; fixed initial variances are not REML estimates. |
| Methylation association | Match CpG/sample matrices to phenotype and covariates; use limma or an appropriate regression with justified phenotype coding and batch/cell-composition adjustment. | Effect size, uncertainty, raw and multiple-testing-adjusted p-values for all tested CpGs. Ordinal metabolizer labels must not silently become equally spaced numbers. |

## Peaks, contacts, and expression

| Task | Inputs and execution | Checks and outputs |
| --- | --- | --- |
| ChIP-seq peaks | Treatment and control alignments, read layout, effective genome size and q-value; use MACS3 in the corresponding paired/single-end mode. | Native peak tables, duplicate policy, signal/control tracks and replicate QC. Record the version; equivalence to older MACS2 results is not established. |
| ATAC differential accessibility | Peak discovery is only the first step. Build a reproducible consensus peak set, count fragments per biological replicate, and fit the requested design with DESeq2/edgeR. | Sample-by-peak counts, effect sizes and FDR with TSS enrichment, library complexity and fragment-length QC. One treatment BAM versus one control BAM is not a replicated differential analysis. |
| Motif enrichment | Foreground intervals/sequences, reference and matched background; use HOMER. | Motif IDs, enrichment statistics and multiple-testing corrections. Match background length, GC and accessibility as required by the hypothesis. |
| Motif sites | Select the actual JASPAR matrix/version and organism, pseudocounts, background and score threshold; scan both strands using Biopython motifs. | Matrix ID, strand, interval and score/p-value. Do not choose the first same-name matrix or confuse relative score with the scanner's raw threshold. |
| Region overlaps | Use BioV intervals with a common assembly and 0-based half-open coordinates. Define whether counting pairs, unique intervals or union-covered bases. | Per-interval overlap table and clearly defined denominators. Multiple overlaps must not inflate coverage beyond the union length. |
| Hi-C interactions / domains | Use cooler to read the correct resolution and balancing weights, then cooltools for expected contacts, insulation or a stated loop method. | Contact/insulation tracks, domain/loop coordinates and filtering parameters. Never silently switch balanced to raw counts; distance-normalize before interpreting contact enrichment. |
| Gene-set enrichment | Fix organism, ID namespace, gene-set release and the measured/tested gene universe. Use GSEApy/Enrichr for overrepresentation or a ranked method for a complete statistic vector. | All tested sets, overlap genes, p-values/FDR and universe. The supported-library list comes from the service, not a bundled static catalog. |
| NMF | Supply a nonnegative gene-by-sample matrix and normalization provenance; use scikit-learn NMF with rank, seed and convergence criteria. | Both factor matrices with feature/sample IDs, residuals and stability across starts/ranks. Reject negative input rather than taking its absolute value. |
| DDR coexpression networks | Use matched expression/mutation matrices, defined DDR genes and NetworkX after testing correlations with FDR control. | Edges, correlations, uncertainty, centrality and mutation frequencies. Correlation plus unequal mutation frequencies does not establish synthetic lethality. |
| Orthology / ARCHS4 | Use Ensembl Compara/BioMart with explicit source/target species and release; use gget/ARCHS4 with the requested species and expression unit. | Preserve one-to-many mappings, unmapped IDs, experiment IDs and units; do not reduce orthology to a first-hit lookup or call every expression value TPM. |

For single-cell QC, integration, transfer and embeddings, use
`single-cell-annotation`. For alignments, haplotypes and phylogeny, use
`analyze-sequences`. Demographic simulation belongs in `model-biological-systems`.

Method documentation: [SuSiE](https://stephenslab.github.io/susieR/),
[GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2),
[CNVkit](https://cnvkit.readthedocs.io/),
[DESeq2](https://bioconductor.org/packages/DESeq2/),
[cooltools](https://cooltools.readthedocs.io/),
[GSEApy](https://gseapy.readthedocs.io/), and
[sommer](https://cran.r-project.org/package=sommer).
Select the method and reference database for the actual task; this skill
does not establish that every analysis has been validated.

Migration source: Biomni commit `400c1f366b96a35ca253e13c9b06c5076af41d65`,
`genomics.py`, `genetics.py`, `cancer_biology.py`, and the ATAC/MWAS functions
in `immunology.py` and `pharmacology.py`. These instructions replace task
decisions and reject the noted scientific shortcuts; they do not preserve old APIs.

Run Python analyses with `biov exec python -- analysis.py` and R analyses with
`biov exec r -- analysis.R`. Invoke native tools such as GATK with
`biov exec gatk -- ...`, retaining their documented arguments. For genomic
prediction, retain the selected sommer method and version; if unavailable,
report that limitation rather than substituting an unvalidated estimator.
Other native tools are `crossmap`, `samtools`, `bcftools`, `snpeff`, `lumpy`
(lumpyexpress), `cnvkit` and `macs3`. The `homer` entry already runs
`annotatePeaks.pl`; select another HOMER program with native Pixi arguments and
the same manifest, for example
`biov pixi run --manifest-path "$BIOV_ENVIRONMENT_MANIFEST" --environment homer -- findMotifsGenome.pl peaks.txt hg38 output/ -size given`.

Read and write files in the working directory or configured data directories.
For a specialized tool without BioV, first obtain input files with BioV in a
`biov run prepare.py` script in the BioV application environment, then pass their paths to the tool.
