# Typed sequence Series

BioV provides three explicit nullable pandas extension dtypes:

```python
import pandas as pd
import biov  # registers BioV extension dtypes and the .seq accessor

dna = pd.Series(["ACGT", None], dtype="biov.dna")
rna = pd.Series(["AUGGCU"], dtype="biov.rna")
protein = pd.Series(["MKWVTF"], dtype="biov.protein")
```

The dtype is required. BioV never guesses sequence kind from an ordinary string Series, and accessing `.seq` on one raises `AttributeError`.

Values are normalized to uppercase. `None` and `pd.NA` remain missing, while an empty string remains a present sequence of length zero. DNA accepts `ACGTRYSWKMBDHVN`; RNA accepts `ACGURYSWKMBDHVN`; protein accepts the 20 canonical amino acids plus `B`, `J`, `O`, `U`, `X`, `Z`, and `*`. Any other symbol or non-string value raises `SequenceValidationError` when the array is created or assigned.

All three types expose nullable integer `Series.seq.length`.

DNA and RNA expose:

- `reverse_complement()`, returning the same sequence dtype;
- `gc_fraction()`, using Biopython's weighted IUPAC ambiguity model and returning 0 for an empty sequence;
- `translate(table=1, to_stop=False)`, returning `biov.protein` and requiring each present nucleotide sequence to contain complete codons.

Protein exposes:

- `molecular_weight()`;
- `isoelectric_point()`;
- `amino_acid_composition()`, a DataFrame containing percentage columns in canonical amino-acid order: `ACDEFGHIKLMNPQRSTVWY`.

Biopython defines protein analyses only for non-empty, unambiguous canonical sequences. BioV still stores extended IUPAC symbols, translated stops, and empty protein strings, but calling one of the three protein analyses on such a value raises `SequenceValidationError` instead of leaking backend-specific exceptions. Missing inputs produce missing outputs and do not raise.
