# Native sequence normalization and reverse complement

This contract defines the first Rust sequence slice: batch normalization and
DNA/RNA reverse complement. It settles these two operations' public Python
representation and error mapping. It does **not** settle the future dataframe
API, migrate scalar `Seq`, or complete the broader SPEC T58 inventory.

## Public Python surface

```python
import biov

biov.normalize_sequences(["acgtryn", None, "", "acgtryn"], kind="dna")
# ["ACGTRYN", None, "", "ACGTRYN"]

biov.reverse_complements(["acgtryn", None, "", "acgtryn"], kind="dna")
# ["NRYACGT", None, "", "NRYACGT"]

biov.reverse_complements(["augcryswkmbdhvn"], kind="rna")
# ["NBDHVKMWSRYGCAU"]
```

Signatures:

```python
normalize_sequences(values: list[str | None], *, kind: str) -> list[str | None]
reverse_complements(values: list[str | None], *, kind: str) -> list[str | None]
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
- Results contain uppercase ASCII strings and `None`, rather than `Seq`,
  `SeqRecord`, pandas scalars or extension arrays. No annotations, indexes,
  descriptions or other metadata are accepted or silently discarded here.
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

## Errors and atomicity

| Invalid input | Exception |
| --- | --- |
| Unknown string kind, including uppercase or padded names | `ValueError` |
| Non-string kind, missing kind or positional kind | `TypeError` |
| A non-list input container | `TypeError` |
| A list element other than string or `None` | `TypeError` |
| A string containing a disallowed symbol | `biov.SequenceValidationError` |
| Reverse complement with protein kind | `TypeError` |

`SequenceValidationError` is a subclass of `ValueError`, also available through
`biov.seq` for the existing pandas surface. Integers, booleans, bytes, `NaN`,
`pd.NA` and sequence objects are not silently converted to strings or missing
values by this list API. Any invalid element rejects the whole batch; a valid
prefix is never returned. Error-message wording is diagnostic, not a stable
serialization format. No exception precedence is promised when several distinct
input requirements are violated simultaneously.

## Relationship to existing interfaces

The existing pandas `biov.dna`, `biov.rna`, `biov.protein` arrays and `.seq`
accessor remain during this slice. Their normalization and reverse-complement
operations delegate to this native computation; this does not add a separate
compatibility layer or retain duplicate scientific algorithms. Pandas still has
its own missing-value and container adaptation. The list API deliberately does
not inherit that adaptation.

For new batch code, replace a pandas construction performed solely to run one of
these two operations with an explicit list and `kind`. Keep any record IDs,
annotations or other columns separately, preserving their association with row
positions. Do not use `str(record)` as an implicit metadata migration.

ASCII-only normalization is an explicit tightening of the old use of Python
`str.upper()`, which could turn certain non-ASCII inputs into accepted symbols.
The scalar `biov.Seq` remains Biopython-backed and is outside this migration.
Translation, weighted GC, protein properties, FASTA/SeqRecord behavior, interval
operations, Polars selection and CLI/MCP transport remain separate work.

## Scientific provenance and acceptance

The biological definitions are checked independently of the new Rust source:

1. [INSDC Feature Table Definition](https://www.insdc.org/submitting-standards/feature-table/),
   version 11.4 (April 2026), §7.4.1 supplies nucleotide ambiguity-set meanings;
   §7.4.3 supplies amino-acid symbols. The nucleotide section cites
   Cornish-Bowden, *Nucleic Acids Research* 13, 3021–3030 (1985).
   Its feature notation uses `t` for RNA uracil; BioV deliberately spells RNA
   uracil `U` in the declared RNA alphabet. These are symbol definitions, not
   copied implementation code.
2. [Biopython 1.88 IUPAC data](https://github.com/biopython/biopython/blob/biopython-188/Bio/Data/IUPACData.py)
   and [sequence implementation](https://github.com/biopython/biopython/blob/biopython-188/Bio/Seq.py)
   are the separately implemented differential reference. `uv.lock` records the
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

These tests are offline after installation. Run against the installed wheel from
outside its source tree, with `PYTHONPATH` unset, as part of the distribution gate;
repeat against a wheel built from an unpacked sdist. Merely passing the Python
source tests is not evidence of either distribution gate. Compilation, fixture
agreement and involution alone are insufficient scientific validation, and no
performance improvement is asserted by this slice. See the
[Rust migration guide](rust-migration.md) for remaining acceptance requirements.
