//! Declared alphabets and IUPAC complementation; see sequence-contract.md.
use std::{collections::BTreeSet, error::Error, fmt, str::FromStr};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Dna,
    Rna,
    Protein,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SequenceError {
    UnknownKind(String),
    InvalidSymbols { kind: Kind, symbols: String },
    NotNucleic,
}

impl fmt::Display for SequenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownKind(kind) => write!(f, "Unknown BioV sequence kind: {kind:?}"),
            Self::InvalidSymbols { kind, symbols } => {
                let label = match kind {
                    Kind::Dna => "DNA",
                    Kind::Rna => "RNA",
                    Kind::Protein => "protein",
                };
                write!(f, "invalid {label} symbol(s): {symbols}")
            }
            Self::NotNucleic => write!(f, "This operation requires DNA or RNA"),
        }
    }
}
impl Error for SequenceError {}

impl FromStr for Kind {
    type Err = SequenceError;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "dna" => Ok(Self::Dna),
            "rna" => Ok(Self::Rna),
            "protein" => Ok(Self::Protein),
            _ => Err(SequenceError::UnknownKind(value.to_owned())),
        }
    }
}
impl Kind {
    fn alphabet(self) -> &'static str {
        match self {
            Self::Dna => "ACGTRYSWKMBDHVN",
            Self::Rna => "ACGURYSWKMBDHVN",
            Self::Protein => "ACDEFGHIKLMNPQRSTVWYBJOUXZ*",
        }
    }
}

/// Normalize ASCII case only. Never transliterate Unicode or infer an alphabet.
pub fn normalize(value: &str, kind: Kind) -> Result<String, SequenceError> {
    let normalized = value.to_ascii_uppercase();
    let invalid: BTreeSet<char> = normalized
        .chars()
        .filter(|symbol| !kind.alphabet().contains(*symbol))
        .collect();
    if !invalid.is_empty() {
        return Err(SequenceError::InvalidSymbols {
            kind,
            symbols: invalid
                .into_iter()
                .map(|c| format!("{c:?}"))
                .collect::<Vec<_>>()
                .join(", "),
        });
    }
    Ok(normalized)
}

/// Reverse and complement validated IUPAC nucleotide sets, preserving ambiguity.
pub fn reverse_complement(value: &str, kind: Kind) -> Result<String, SequenceError> {
    if kind == Kind::Protein {
        return Err(SequenceError::NotNucleic);
    }
    let normalized = normalize(value, kind)?;
    Ok(normalized
        .bytes()
        .rev()
        .map(|symbol| match symbol {
            b'A' if kind == Kind::Rna => 'U',
            b'A' => 'T',
            b'T' | b'U' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            b'R' => 'Y',
            b'Y' => 'R',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'V' => 'B',
            b'D' => 'H',
            b'H' => 'D',
            b'S' => 'S',
            b'W' => 'W',
            b'N' => 'N',
            _ => unreachable!("normalize accepts only the declared nucleotide alphabet"),
        })
        .collect())
}

/// One batch call with stable order, duplicates and distinct null/empty values.
pub fn normalize_batch(
    values: &[Option<String>],
    kind: Kind,
) -> Result<Vec<Option<String>>, SequenceError> {
    values
        .iter()
        .map(|value| value.as_deref().map(|s| normalize(s, kind)).transpose())
        .collect()
}

pub fn reverse_complement_batch(
    values: &[Option<String>],
    kind: Kind,
) -> Result<Vec<Option<String>>, SequenceError> {
    if kind == Kind::Protein {
        return Err(SequenceError::NotNucleic);
    }
    values
        .iter()
        .map(|value| {
            value
                .as_deref()
                .map(|s| reverse_complement(s, kind))
                .transpose()
        })
        .collect()
}

/// Validated biological-symbol count, not UTF-8 bytes of unchecked input.
pub fn length(value: &str, kind: Kind) -> Result<usize, SequenceError> {
    Ok(normalize(value, kind)?.len())
}

/// Mean GC probability with equal weight for the bases represented by each code.
///
/// All symbols contribute to the denominator. Empty sequences return zero by
/// BioV's explicit convention. Sixths make the IUPAC weights exact until the
/// final f64 division; a u128 sum cannot overflow for an addressable Rust string.
pub fn weighted_gc_fraction(value: &str, kind: Kind) -> Result<f64, SequenceError> {
    Ok(nucleotide_gc_counts(value, kind)?.weighted_gc_fraction())
}

/// Exact additive nucleotide GC counts. The weighted numerator uses sixths,
/// preserving exact IUPAC probabilities when independently read chunks are joined.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NucleotideGcCounts {
    pub length: u64,
    pub canonical_base_count: u64,
    pub gc_base_count: u64,
    pub weighted_gc_sixths: u128,
}

impl NucleotideGcCounts {
    /// Literal G/C among only canonical A/C/G/T (or A/C/G/U for RNA).
    /// An all-ambiguous or empty sequence has no canonical fraction.
    pub fn gc_fraction(self) -> Option<f64> {
        (self.canonical_base_count != 0)
            .then(|| self.gc_base_count as f64 / self.canonical_base_count as f64)
    }

    /// Equal-weight IUPAC GC probability among all symbols; empty is zero.
    pub fn weighted_gc_fraction(self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            self.weighted_gc_sixths as f64 / (6.0 * self.length as f64)
        }
    }
}

/// Validate the declared alphabet and compute reusable exact GC numerators.
pub fn nucleotide_gc_counts(value: &str, kind: Kind) -> Result<NucleotideGcCounts, SequenceError> {
    if kind == Kind::Protein {
        return Err(SequenceError::NotNucleic);
    }
    let normalized = normalize(value, kind)?;
    let mut counts = NucleotideGcCounts {
        length: normalized.len() as u64,
        ..NucleotideGcCounts::default()
    };
    for symbol in normalized.bytes() {
        if matches!(symbol, b'A' | b'C' | b'G' | b'T' | b'U') {
            counts.canonical_base_count += 1;
        }
        if matches!(symbol, b'C' | b'G') {
            counts.gc_base_count += 1;
        }
        counts.weighted_gc_sixths += match symbol {
            b'G' | b'C' | b'S' => 6_u128,
            b'A' | b'T' | b'U' | b'W' => 0,
            b'R' | b'Y' | b'K' | b'M' | b'N' => 3,
            b'B' | b'V' => 4,
            b'D' | b'H' => 2,
            _ => unreachable!("normalize accepts only the declared nucleotide alphabet"),
        };
    }
    Ok(counts)
}

/// Compute validated lengths in stable batch order with nulls intact.
pub fn length_batch(
    values: &[Option<String>],
    kind: Kind,
) -> Result<Vec<Option<usize>>, SequenceError> {
    values
        .iter()
        .map(|value| value.as_deref().map(|s| length(s, kind)).transpose())
        .collect()
}

/// Compute weighted GC fractions in stable batch order with nulls intact.
pub fn weighted_gc_fraction_batch(
    values: &[Option<String>],
    kind: Kind,
) -> Result<Vec<Option<f64>>, SequenceError> {
    if kind == Kind::Protein {
        return Err(SequenceError::NotNucleic);
    }
    values
        .iter()
        .map(|value| {
            value
                .as_deref()
                .map(|s| weighted_gc_fraction(s, kind))
                .transpose()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn full_alphabets() {
        assert_eq!(
            reverse_complement("ACGTRYSWKMBDHVN", Kind::Dna).unwrap(),
            "NBDHVKMWSRYACGT"
        );
        assert_eq!(
            reverse_complement("ACGURYSWKMBDHVN", Kind::Rna).unwrap(),
            "NBDHVKMWSRYACGU"
        );
        assert_eq!(
            normalize("acdefghiklmnpqrstvwybjouxz*", Kind::Protein).unwrap(),
            "ACDEFGHIKLMNPQRSTVWYBJOUXZ*"
        );
    }
    #[test]
    fn rejection_and_empty() {
        for value in ["U", "-", " ", "ß", "K", "é", "A\n", "\0", "1"] {
            assert!(normalize(value, Kind::Dna).is_err(), "{value:?}");
        }
        assert!(normalize("T", Kind::Rna).is_err());
        assert!(reverse_complement_batch(&[], Kind::Protein).is_err());
        assert_eq!(reverse_complement("", Kind::Dna).unwrap(), "");
        assert!("DNA".parse::<Kind>().is_err());
    }
    #[test]
    fn batches_and_involution() {
        let values = [
            Some("acgtryswkmbdhvn".into()),
            None,
            Some("".into()),
            Some("aa".into()),
            Some("aa".into()),
        ];
        for kind in [Kind::Dna, Kind::Rna] {
            let values: Vec<_> = values
                .iter()
                .map(|s| {
                    s.as_ref().map(|s: &String| {
                        if kind == Kind::Rna {
                            s.replace('t', "u")
                        } else {
                            s.clone()
                        }
                    })
                })
                .collect();
            let reverse = reverse_complement_batch(&values, kind).unwrap();
            assert_eq!(
                reverse_complement_batch(&reverse, kind).unwrap(),
                normalize_batch(&values, kind).unwrap()
            );
            assert_eq!(reverse[1], None);
            assert_eq!(reverse[2], Some(String::new()));
            assert_eq!(reverse[3], reverse[4]);
        }
    }
    #[test]
    fn sequence_metrics() {
        // Hand-derived GC probabilities for each IUPAC set in alphabet order.
        let expected_sixths = [0, 6, 6, 0, 3, 3, 6, 0, 3, 3, 4, 2, 2, 4, 3];
        for (symbol, weight) in "ACGTRYSWKMBDHVN".chars().zip(expected_sixths) {
            assert_eq!(
                weighted_gc_fraction(&symbol.to_string(), Kind::Dna).unwrap(),
                weight as f64 / 6.0
            );
        }
        assert_eq!(weighted_gc_fraction("GCN", Kind::Dna).unwrap(), 5.0 / 6.0);
        assert_eq!(weighted_gc_fraction("GDVV", Kind::Dna).unwrap(), 2.0 / 3.0);
        assert_eq!(weighted_gc_fraction("au", Kind::Rna).unwrap(), 0.0);
        assert_eq!(weighted_gc_fraction("", Kind::Dna).unwrap(), 0.0);
        assert_eq!(length("m*x", Kind::Protein).unwrap(), 3);
        assert!(length("U", Kind::Dna).is_err());
        assert!(length("ß", Kind::Protein).is_err());
        assert!(weighted_gc_fraction("T", Kind::Rna).is_err());
        assert!(weighted_gc_fraction_batch(&[], Kind::Protein).is_err());
        let values = [None, Some("".into()), Some("gcN".into())];
        assert_eq!(
            length_batch(&values, Kind::Dna).unwrap(),
            vec![None, Some(0), Some(3)]
        );
        assert_eq!(
            weighted_gc_fraction_batch(&values, Kind::Dna).unwrap(),
            vec![None, Some(0.0), Some(5.0 / 6.0)]
        );
    }

    #[test]
    fn exact_gc_counts_preserve_ambiguity_denominators_and_chunk_additivity() {
        let all = nucleotide_gc_counts("acgtryswkmbdhvnGCN", Kind::Dna).unwrap();
        assert_eq!(all.length, 18);
        assert_eq!(all.canonical_base_count, 6);
        assert_eq!(all.gc_base_count, 4);
        assert_eq!(all.weighted_gc_sixths, 60);
        assert_eq!(all.gc_fraction(), Some(2.0 / 3.0));
        assert_eq!(all.weighted_gc_fraction(), 5.0 / 9.0);
        let a = nucleotide_gc_counts("acgtrysw", Kind::Dna).unwrap();
        let b = nucleotide_gc_counts("kmbdhvnGCN", Kind::Dna).unwrap();
        assert_eq!(
            a.weighted_gc_sixths + b.weighted_gc_sixths,
            all.weighted_gc_sixths
        );
        assert_eq!(
            a.canonical_base_count + b.canonical_base_count,
            all.canonical_base_count
        );
        let ambiguous = nucleotide_gc_counts("nSWbDhV", Kind::Dna).unwrap();
        assert_eq!(ambiguous.gc_fraction(), None);
        assert_eq!(ambiguous.weighted_gc_fraction(), 0.5);
        assert_eq!(
            nucleotide_gc_counts("", Kind::Dna).unwrap().gc_fraction(),
            None
        );
        assert_eq!(
            nucleotide_gc_counts("uGc", Kind::Rna)
                .unwrap()
                .gc_fraction(),
            Some(2.0 / 3.0)
        );
        assert!(nucleotide_gc_counts("uGc", Kind::Dna).is_err());
        assert!(nucleotide_gc_counts("G", Kind::Protein).is_err());
    }
}
