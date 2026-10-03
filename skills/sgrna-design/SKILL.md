---
name: sgrna-design
description: Design CRISPR sgRNAs with a three-tier approach — validated Addgene sequences first, pre-computed CRISPick datasets second, de novo design rules last. Use whenever asked to find, select, or design guide RNAs for knockout, CRISPRa, or CRISPRi experiments.
license: CC BY 4.0
---

# sgRNA Design

Follow the three-tier approach. Always start at tier 1 and escalate only when
the previous tier yields nothing.

1. **Validated sequences first.** Search `references/addgene_grna_sequences.csv`
   (300+ experimentally validated guides) by `Target_Gene`, `Target_Species`,
   and `Application`. Also run a literature/web search before giving up — many
   validated sgRNAs are published but not in Addgene.
2. **CRISPick pre-computed designs.** Pick the right dataset from
   `references/CRISPick_download_links.txt` by taxonomy ID, genome build, Cas
   enzyme, and application; download, filter for the gene, and rank by
   `Combined Rank`. AsCas12a and enAsCas12a are different enzymes — never mix
   their datasets.
3. **De novo design rules.** Only when the gene or organism is not covered:
   correct guide length and PAM for the enzyme, 40-60% GC, no TTTT or >4
   homopolymer runs, early exons for knockout, TSS windows for CRISPRa/i.

Always recommend testing 3-4 sgRNAs per gene, and cite the PubMed ID attached
to any validated sequence used.

Read `references/sgrna-design-guide.md` for the complete workflow, file
formats, ranking columns, filtering recipes, and worked examples.
