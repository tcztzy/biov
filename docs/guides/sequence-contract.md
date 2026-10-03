# Native sequence normalization, reverse complement and metrics

This contract defines the native Rust sequence slices: batch normalization,
DNA/RNA reverse complement, validated sequence length and weighted-IUPAC GC
fraction. It settles these operations' public Python representation and error
mapping. It does **not** settle the future dataframe API, migrate scalar `Seq`,
or complete the broader SPEC T58 inventory.

## Public Python surface

```python
import biov

biov.normalize_sequences(["acgtryn", None, "", "acgtryn"], kind="dna")
# ["ACGTRYN", None, "", "ACGTRYN"]

biov.reverse_complements(["acgtryn", None, "", "acgtryn"], kind="dna")
# ["NRYACGT", None, "", "NRYACGT"]

biov.reverse_complements(["augcryswkmbdhvn"], kind="rna")
# ["NBDHVKMWSRYGCAU"]

biov.sequence_lengths(["acgtryn", None, "", "acgtryn"], kind="dna")
# [7, None, 0, 7]

biov.sequence_lengths(["M*", None, "", "BXZJUO*"], kind="protein")
# [2, None, 0, 7]

biov.weighted_gc_fractions(["GCN", None, "", "NN", "GDVV"], kind="dna")
# [0.8333333333333334, None, 0.0, 0.5, 0.6666666666666666]
```

Signatures:

```python
normalize_sequences(values: list[str | None], *, kind: str) -> list[str | None]
reverse_complements(values: list[str | None], *, kind: str) -> list[str | None]
sequence_lengths(values: list[str | None], *, kind: str) -> list[int | None]
weighted_gc_fractions(values: list[str | None], *, kind: str) -> list[float | None]
```

- `kind` is required and keyword-only. Its accepted values are the exact strings
  `"dna"`, `"rna"` and `"protein"`. There is no alphabet inference.
- The input is a Python list of strings or `None`. Iterables, pandas objects,
  NumPy arrays, scalar strings and bytes are not this API's input representation.
  Convert a container deliberately before calling; do not stringify its elements.
- Each call returns a new list. Original row order, duplicates, list length and
  null positions are preserved. Inputs are never mutated, including on failure.
- `None` is the sole missing value. An empty string is a present, zero-length
  sequence. Empty lists and all-null lists are valid for supported operations.
- Normalization and reverse complement return uppercase ASCII strings and
  `None`; lengths return Python integers and `None`; GC fractions return Python
  floats and `None`. Results never contain `Seq`, `SeqRecord`, pandas scalars or
  extension arrays. No annotations, indexes, descriptions or other metadata are
  accepted or silently discarded here.
- Validation and computation cross into the native extension once per batch.
  There is no Python algorithm fallback when the extension is unavailable.

## Alphabet and normalization

| Kind | Complete accepted uppercase alphabet |
| --- | --- |
| DNA | `ACGTRYSWKMBDHVN` |
| RNA | `ACGURYSWKMBDHVN` |
| Protein | `ACDEFGHIKLMNPQRSTVWYBJOUXZ*` |

Only the ASCII lowercase counterparts of these letters are additionally accepted;
normalization replaces them with uppercase. `*` is retained as a protein stop
marker. Extended protein symbols are storable here; accepting them does not imply
that every protein-property algorithm accepts them.

Reject non-ASCII characters, whitespace (including leading/trailing spaces,
newlines and tabs), alignment gaps (`-` or `.`), digits and all other characters
outside the declared alphabet, including lone surrogate code points in Python
strings. DNA rejects `U/u`; RNA rejects `T/t`. Neither is
silently transcribed. Nucleotide `X/x` is rejected even though some reference
libraries accept it. No whitespace is stripped, gaps removed, residues replaced
or Unicode case-folding applied. In particular, `ß`, `ſ` and `ı` must not become
valid multi-letter or ASCII sequences through Python Unicode uppercasing.

## Reverse-complement definition

Input and output are both written 5′ to 3′. Reverse the symbol order, then
complement each symbol. For an ambiguity symbol representing a set of bases,
complement every member and encode the resulting set using the matching IUPAC
symbol. DNA uses A↔T and C↔G; RNA uses A↔U and C↔G.

| Input symbol | DNA complement | RNA complement |
| --- | --- | --- |
| A | T | U |
| C | G | G |
| G | C | C |
| T | A | invalid |
| U | invalid | A |
| R | Y | Y |
| Y | R | R |
| S | S | S |
| W | W | W |
| K | M | M |
| M | K | K |
| B | V | V |
| D | H | H |
| H | D | D |
| V | B | B |
| N | N | N |

Complete, hand-derived examples used as acceptance fixtures:

| Kind | Input | Complete output |
| --- | --- | --- |
| DNA | `ACGTRYSWKMBDHVN` | `NBDHVKMWSRYACGT` |
| RNA | `ACGURYSWKMBDHVN` | `NBDHVKMWSRYACGU` |
| DNA | `aCGTryn` | `NRYACGT` |
| RNA | `aCGUryn` | `NRYACGU` |
| DNA | `AGTC` | `GACT` |
| RNA | `AGUC` | `GACU` |

`reverse_complements` supports only DNA and RNA. Declaring protein raises
`TypeError` even for `[]`, `[None]` or `[""]`; no sequence content can make this
operation meaningful for a protein batch.

## Validated sequence-length definition

`sequence_lengths` supports DNA, RNA and protein. It counts the number of
symbols after validating the complete sequence against its declared alphabet.
ASCII lowercase letters count exactly like their uppercase counterparts. Every
accepted symbol counts once, including nucleotide ambiguities, extended protein
symbols and each protein stop marker `*`. A stop marker does not terminate the
count. This is a symbol count, not a codon, translated-residue or molecular-length
calculation. Invalid symbols are rejected even though their string length could
otherwise be computed.

The empty string returns integer `0`; `None` returns `None`. Complete examples:

| Kind | Input | Complete output |
| --- | --- | --- |
| DNA | `["acgtryswkmbdhvn", "", None, "GCN"]` | `[15, 0, None, 3]` |
| RNA | `["acguryswkmbdhvn", "", None, "GCN"]` | `[15, 0, None, 3]` |
| Protein | `["M*", "***", "B*Z*", "BXZJUO*", "", None]` | `[2, 3, 4, 7, 0, None]` |

## Weighted-IUPAC GC definition

`weighted_gc_fractions` supports only DNA and RNA. Each ambiguity symbol means
an equally weighted set of possible canonical bases. Its GC contribution is the
number of `G` or `C` members divided by the total number of members in that set.
For a nonempty sequence, sum those contributions and divide by the number of
sequence symbols. **Every accepted position remains in the denominator**,
including `A`, `T`/`U`, `W` and all ambiguous positions. This is a fraction from
`0.0` to `1.0`, not a percentage or an estimate using observed base frequencies.

| Symbol(s) | Exact GC contribution per position |
| --- | --- |
| `A`, `T` (DNA), `U` (RNA), `W` | `0` |
| `C`, `G`, `S` | `1` |
| `R`, `Y`, `K`, `M`, `N` | `1/2` |
| `B`, `V` | `2/3` |
| `D`, `H` | `1/3` |

The result for each present sequence is a finite Python float, subject to ordinary
floating-point rounding. The empty string returns `0.0`; `None` returns `None`.
Ambiguous positions are neither removed nor treated as entirely non-GC.
Complete hand-derived examples, valid for both DNA and RNA:

| Input | Exact result (returned as float) |
| --- | --- |
| `GCN` | `(1 + 1 + 1/2) / 3 = 5/6` |
| `GDVV` | `(1 + 1/3 + 2/3 + 2/3) / 4 = 2/3` |
| `AWN` | `(0 + 0 + 1/2) / 3 = 1/6` |
| `SW` | `(1 + 0) / 2 = 1/2` |
| `NNNN` | `1/2` |
| `""` | `0.0` |

Reverse complement preserves weighted GC because canonical complement pairing
preserves membership in the `G`/`C` set. Declaring protein raises `TypeError`,
even for `[]`, `[None]`, an all-null batch or `[""]`. No gap, `X`, mixed DNA/RNA,
or other out-of-alphabet exception is introduced for GC analysis.

## Errors and atomicity

| Invalid input | Exception |
| --- | --- |
| Unknown string kind, including uppercase or padded names | `ValueError` |
| Non-string kind, missing kind or positional kind | `TypeError` |
| A non-list input container | `TypeError` |
| A list element other than string or `None` | `TypeError` |
| A string containing a disallowed symbol | `biov.SequenceValidationError` |
| Reverse complement or weighted GC with protein kind | `TypeError` |

`SequenceValidationError` is a subclass of `ValueError`, also available through
`biov.seq` for the existing pandas surface. Integers, booleans, bytes, `NaN`,
`pd.NA` and sequence objects are not silently converted to strings or missing
values by this list API. Any invalid element rejects the whole batch; a valid
prefix is never returned. Error-message wording is diagnostic, not a stable
serialization format. No exception precedence is promised when several distinct
input requirements are violated simultaneously.

## Relationship to existing interfaces

The existing pandas `biov.dna`, `biov.rna`, `biov.protein` arrays and `.seq`
accessor remain during these slices. Normalization, reverse complement,
`.seq.length` and `.seq.gc_fraction()` delegate to native computation once per
batch; this does not add a separate compatibility layer or retain duplicate
scientific algorithms. The metric adapters preserve the input Series index and
name and return nullable `Int64` length or `Float64` GC Series. Pandas still has
its own missing-value and container adaptation. The list API deliberately does
not inherit that adaptation.

For new batch code, replace a pandas construction performed solely to run one of
these operations with an explicit list and `kind`. Keep any record IDs,
annotations or other columns separately, preserving their association with row
positions. Do not use `str(record)` as an implicit metadata migration.

ASCII-only normalization is an explicit tightening of the old use of Python
`str.upper()`, which could turn certain non-ASCII inputs into accepted symbols.
The scalar `biov.Seq` remains Biopython-backed and is outside this migration.
Translation, protein properties, FASTA/SeqRecord behavior, interval operations,
Polars selection and CLI/MCP transport remain separate work.

## Scientific provenance and acceptance

The biological definitions are checked independently of the new Rust source:

1. [INSDC Feature Table Definition](https://www.insdc.org/submitting-standards/feature-table/),
   version 11.4 (April 2026), §7.4.1 supplies nucleotide ambiguity-set meanings;
   §7.4.3 supplies amino-acid symbols. The nucleotide section cites
   Cornish-Bowden, *Nucleic Acids Research* 13, 3021–3030 (1985).
   Its feature notation uses `t` for RNA uracil; BioV deliberately spells RNA
   uracil `U` in the declared RNA alphabet. These are symbol definitions, not
   copied implementation code.
2. [Biopython 1.88 IUPAC data](https://github.com/biopython/biopython/blob/biopython-188/Bio/Data/IUPACData.py),
   [sequence implementation](https://github.com/biopython/biopython/blob/biopython-188/Bio/Seq.py)
   and [weighted GC implementation](https://github.com/biopython/biopython/blob/biopython-188/Bio/SeqUtils/__init__.py)
   are the separately implemented differential references. `uv.lock` records the
   exact 1.88 distribution and hashes. Tests pin the reference version so a
   dependency update cannot silently redefine the comparison baseline. BioV
   imports the installed reference for tests; no Biopython implementation source
   is vendored into this slice. Biopython is distributed under its Biopython
   License Agreement or BSD 3-Clause License.

`tests/test_native_sequence.py` retains complete expected outputs and checks:

- Every nucleotide symbol independently, both DNA/RNA full-alphabet fixtures,
  complete protein normalization and asymmetric examples that detect failure to
  reverse the sequence
- An independent set-based oracle, using published ambiguity meanings and only
  the four canonical-base pairings rather than copying Rust's complement table
- Exhaustive nucleotide words through length three, plus deterministic seeded
  mixed-case and longer batches, compared in full against Biopython 1.88
- Involution, length/null/order/duplicate preservation and batch partitioning
- Every ASCII code point, invalid Unicode (including surrogate strings), Python
  types and kinds, including empty-batch
  wrong-kind rejection and failure without input mutation
- A real compiled `_native` module, and execution without calling Biopython's
  reverse-complement methods

`tests/test_native_metrics.py` adds independent acceptance checks for metrics:

- All nucleotide symbols and protein length symbols, complete hand-derived GC
  outputs and repeated protein stop markers
- An independent rational-number GC oracle derived from IUPAC base sets, rather
  than a copied numeric weight table or implementation-generated fixtures
- Full comparisons for exhaustive nucleotide words through length three, seeded
  mixed-case and long inputs, and repeated third-weight ambiguity symbols,
  against both the rational oracle and pinned Biopython 1.88
- Finite fraction results within `[0, 1]`, reverse-complement invariance,
  batch partitioning, complete row/null preservation and input immutability
- Every ASCII code point, Unicode including surrogates, strict list/type/kind
  validation and unsupported protein GC on empty or all-null batches
- Compiled native public callables, one native call per pandas metric batch,
  nullable dtype/index/name preservation and no Biopython GC fallback execution

Integer length outputs are compared exactly. GC differential and rational-oracle
checks allow `1e-12` absolute and relative tolerance for floating-point rounding;
they compare every result and do not accept only preview prefixes or aggregates.

These tests are offline after installation. Run against the installed wheel from
outside its source tree, with `PYTHONPATH` unset, as part of the distribution gate;
repeat against a wheel built from an unpacked sdist. Merely passing the Python
source tests is not evidence of either distribution gate. Compilation, fixture
agreement and involution alone are insufficient scientific validation, and no
performance improvement is asserted by this slice. See the
[Rust migration guide](rust-migration.md) for remaining acceptance requirements.
