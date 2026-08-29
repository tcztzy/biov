# Genomic ranges

`BioDataFrame` owns its public range behavior. RuRanges is an internal Rust/NumPy kernel, not a public container or compatibility layer. Range methods accept another `BioDataFrame`; BioV does not accept or expose PyRanges objects.

## Coordinates and validation

All intervals use 0-based, end-exclusive coordinates: `[start, end)`. Every row must have a non-null string sequence ID and integer coordinates satisfying `0 <= start < end`. Therefore `[0, 10)` and `[10, 20)` touch but do not overlap.

Column names default to `seqid`, `start`, `end`, and `strand` and can be changed with the corresponding `*_col` arguments. If both operands contain the selected strand column, values must be `+` or `-` and operations group by the exact `(seqid, strand)` pair. If neither contains it, operations group by `seqid`. A strand column on only one operand is an error. Pass `strand_col=None` to ignore strand explicitly.

Results use a fresh `RangeIndex`. They preserve self-row input order. Pair-expanding operations order matches by self row and then other-row input order; subtraction orders fragments by self row and ascending coordinate. Duplicate rows remain distinct inputs.

## Operations

`overlap` returns matching rows from the caller, once per input row:

- `first` and `last` select overlap membership and produce the same self-row schema;
- `containment` requires the self interval to contain an interval in `other`;
- `member` requires the self interval to be contained by an interval in `other`.

`intersect` returns one caller-metadata row per overlapping pair, with `start` and `end` clipped to the shared half-open interval. Repeated intervals and repeated matches intentionally produce repeated output rows.

`subtract_ranges` removes the union of matching `other` intervals and returns every non-empty residual fragment with the caller's metadata.

`nearest` appends the selected row from `other`. Colliding other-column names receive `suffix` (default `_b`), the other sequence-ID column is omitted, and `Distance` is added. Overlap has distance 0. Otherwise distance is the uncovered gap plus one, so adjacent half-open intervals have distance 1. At most one result is returned per caller row; a tie selects the lowest other input row, and a row with no candidate in its group is omitted.

Nearest directions are explicit:

| `how` | Meaning |
|---|---|
| `next` | genomic right, independent of strand |
| `previous` | genomic left, independent of strand |
| `upstream` | 5-prime direction: left on `+`, right on `-` |
| `downstream` | 3-prime direction: right on `+`, left on `-` |

`upstream` and `downstream` require valid strand columns on both operands. Overlapping candidates remain eligible in every direction and win with distance 0.

Empty inputs never enter the Rust kernel. Selection/intersection return an empty caller schema, subtraction by an empty operand returns all caller rows, and nearest returns an empty joined schema when either side is empty.
