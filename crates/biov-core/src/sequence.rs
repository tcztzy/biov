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
            Self::NotNucleic => write!(f, "Reverse complement requires DNA or RNA"),
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
}
