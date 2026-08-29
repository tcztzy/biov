# SPEC

## §G GOAL
Replace PyRanges 0.x with BioV-owned RuRanges interval adaptation and add explicit typed pandas DNA/RNA/protein `.seq` APIs.

## §C CONSTRAINTS
- Python `>=3.12`; RuRanges only interval kernel; Biopython only sequence-algorithm backend
- interval semantics ≠ sequence semantics; ⊥ new DataFrame subclass or large framework
- ⊥ PyRanges objects, conversion paths, compatibility aliases, or legacy interval dependencies
- public behavior ! documented & acceptance-tested

## §I INTERFACES
- api: `BioDataFrame.overlap(other, how, seqid_col, start_col, end_col, strand_col)` → selected self rows
- api: `BioDataFrame.intersect(other, seqid_col, start_col, end_col, strand_col)` → clipped self rows per overlap pair
- api: `BioDataFrame.subtract_ranges(other, seqid_col, start_col, end_col, strand_col)` → residual self fragments
- api: `BioDataFrame.nearest(other, seqid_col, start_col, end_col, strand_col, suffix, how)` → self + nearest other columns + `Distance`
- dtype: `biov.dna` | `biov.rna` | `biov.protein` → validated nullable uppercase sequence storage
- accessor: `Series.seq.length` → nullable integer Series
- accessor: DNA/RNA `Series.seq.reverse_complement()` | `gc_fraction()` | `translate(table=1, to_stop=False)`
- accessor: protein `Series.seq.molecular_weight()` | `isoelectric_point()` → nullable float Series; `amino_acid_composition()` → 20-column percentage DataFrame

## §R RESEARCH
id|topic|finding|src
R1|RuRanges surface|stateless `ruranges.numpy` functions accept/return NumPy arrays; groups integer-coded; strand boolean only where required|https://github.com/pyranges/ruranges_py
R2|RuRanges kernels|`overlaps`, `nearest`, `subtract` return source indices; nearest physical directions = `forward`/`backward` & overlap distance = 0|https://raw.githubusercontent.com/pyranges/ruranges_py/master/ruranges/numpy.py
R3|pandas extension|custom dtype + 1-D ExtensionArray preserve semantic type; accessor init ! reject wrong dtype with `AttributeError`|https://pandas.pydata.org/docs/development/extending.html
R4|Biopython sequence math|weighted GC defines ambiguous IUPAC handling & empty GC = 0; molecular weight requires unambiguous residues|https://biopython.org/docs/latest/api/Bio.SeqUtils.html

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

## §T TASKS
id|status|task|cites
T1|x|write contracts & failing acceptance tests|I.*,V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13
T2|x|replace interval adapter with RuRanges NumPy kernels|I.overlap,I.intersect,I.subtract_ranges,I.nearest,V1,V2,V3,V4,V5,V6,V7,V8,V9,V14
T3|x|add typed sequence EA/dtypes/accessor|I.dtype,I.accessor,V10,V11,V12,V13,V14
T4|x|update public docs, exports, dependencies & lock|V9,V15
T5|x|run full verification matrix|V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15

## §B BUGS
id|date|cause|fix
