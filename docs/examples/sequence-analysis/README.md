# Complete CDS translation and protein properties

This example checks annotated CDS translation, then calculates protein sequence
properties. It uses Biopython directly and has no BioV dependency inside its
scientific environment. These calculations do not establish cold tolerance,
protein activity, or experimentally measured molecular weight or pI.

## Execution contract

Both scripts receive two arguments: an `inputs.json` file mapping input names to
absolute local file paths, and a `parameters.json` file. Run each step in a new
directory. Each writes `checks.json`, recording checks actually performed; a
failed check has `passed: false` and the process exits unsuccessfully.

- `extract_cds.py`: input name `genbank`; parameter `accessions` is the explicit
  list of records to analyze. Produces `proteins.fasta` and `cds.csv`.
- `protein_properties.py`: input name `proteins`, pointing to the complete first
  step FASTA; parameters `{}`. Produces `properties.csv`.

For the bundled example select
`["X55053.1", "X62281.1", "M81224.1", "L31939.1", "AF297471.1"]`.
All five proteins must reach the second step, independently of any preview limit.
The complete results retain their input order. A preview of two rows therefore
omits three records without changing either analysis.

Install the native Pixi environment with:

```sh
pixi install --manifest-path docs/examples/sequence-analysis/pyproject.toml --environment python --locked
```

The manifest requires Pixi 0.81.0 and locks Python 3.12 and Biopython 1.88 for
`osx-arm64` and `linux-64`. Select this manifest when executing the scripts through
BioV; its `python` environment runs the ordinary Python executable.

## Scientific scope and checks

The first step requires exactly one complete CDS for each selected record,
standard genetic code 1, `codon_start=1`, and no translation exceptions. It uses
`SeqFeature.translate(..., cds=True)` to handle joined features and strand, check
the complete-CDS rules, and compare the result with the annotated translation.
It preserves accession versions and protein IDs. CSV locations are Biopython's
0-based, end-exclusive intervals; GenBank's original feature coordinates are
1-based and inclusive. Joined locations retain every part. Each location belongs
to its own source accession, not a shared genome assembly.

The second step requires unique identifiers and nonempty sequences containing
only the canonical 20 amino acids. `ProteinAnalysis(monoisotopic=False)` computes
average molecular weight in daltons and theoretical isoelectric point. Modified
residues, ambiguous residues, and stop symbols are rejected, not removed.

The five protein lengths are 66, 67, 65, 65, and 65 amino acids. Their translations
match the supplied GenBank annotations. The sixth fixture record, `AJ237582.1`,
has a partial CDS and is intentionally not selected. Selecting it must fail the
complete-CDS requirement. The fixture includes two joined complete CDSs, which
exercise spliced extraction. It is a small validation example, not a general
annotation pipeline; alternative codes and partial CDS analyses need an explicit
different method.

Run the independent scientific checks in the BioV development environment:

```sh
python -m pytest tests/test_analysis_example.py
```

## Fixture source and license

`tests/data/cor6_6.gb` is copied unchanged from the
[Biopython 1.88 GenBank test fixture](https://raw.githubusercontent.com/biopython/biopython/biopython-188/Tests/GenBank/cor6_6.gb),
tag `biopython-188`. The original records retain their accession versions,
organisms, annotations, and literature references.

- Size: 14,967 bytes.
- SHA-256: `01b4e193b71344752a96e2b37e886118a33c07b16c3684741c31eb3efeaf60b8`.
- Distribution terms: [upstream LICENSE.rst at the same tag](https://github.com/biopython/biopython/blob/biopython-188/LICENSE.rst).

Biopython's upstream license states that files without a different individual
license header use the Biopython License Agreement. Its notice is reproduced
below for the copied fixture. Biopython contributors retain their notices; no
endorsement is implied.

### Biopython License Agreement

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.
