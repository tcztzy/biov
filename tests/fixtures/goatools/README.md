# Synthetic GO enrichment acceptance case

These are artificial IDs and a three-term local ontology, not a GO database
snapshot or biological interpretation. No remote annotations or ontology are
read. `annotations.id2gos` assigns genes 1–4 to alpha and genes 5–10 to beta;
all genes are in the population, and genes 1–4 form the study. GOATOOLS propagates
these counts to the biological-process root via `is_a`.

`expected.tsv` is the complete native CLI output from the clean Linux-64
installation documented in `docs/guides/environments.md`. It preserves all three
tested terms, full floating-point p-values, counts, and study members.

Settings: GOATOOLS 1.6.5, Python 3.12.14, SciPy 1.18.1, statsmodels 0.14.6;
BP only, two-sided `fisher_scipy_stats`, alpha 0.05, pval 1 (retain every result),
Bonferroni and Benjamini–Hochberg, default count propagation, `id2gos` format.

An independent exact integer hypergeometric enumeration in
`scripts/validate_goatools.py` checks the CLI output, rather than using GOATOOLS
or SciPy to manufacture its expectations. Alpha has the table [[4,0],[0,6]],
beta [[0,4],[6,0]], root [[4,0],[6,0]]. The leaf p-values are 1/210,
Bonferroni-adjusted 1/70 (three tests), and BH-adjusted 1/140 (two tied smallest
p-values among three tests); all root p-values are 1.
