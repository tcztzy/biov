//! Offline biological identifier syntax for a deliberately curated subset.
//!
//! This crate validates RefSeq GCF assembly references and UniProtKB accession
//! references; it does not resolve them, check that records exist, or establish
//! that a UniProt accession is primary rather than secondary. It is not a port
//! of Python BioV's complete identifiers.org registry or prompt scanner.
//!
//! A biological identifier names an external biological record. An opaque
//! BioV dataset handle names a locally managed result and must remain a separate
//! type: `dataset://...`, paths, and arbitrary handle strings are not identifiers.
//!
//! ```
//! use biov_identifiers::{IdentifierRef, Namespace};
//! let assembly = IdentifierRef::parse("refseq.gcf://GCF_000001405.40")?;
//! assert_eq!(assembly.namespace(), Namespace::RefSeqGcf);
//! assert_eq!(assembly.base_accession(), "GCF_000001405");
//! assert_eq!(assembly.version(), Some("40"));
//! assert_eq!(assembly.compact_id(), "refseq.gcf:GCF_000001405.40");
//! # Ok::<(), biov_identifiers::IdentifierError>(())
//! ```

use serde::{Deserialize, Serialize};
use std::{error::Error, fmt, str::FromStr};

/// Namespaces whose syntax is supported locally. No automatic registry lookup.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Namespace {
    /// NCBI RefSeq assembly; distinct from GenBank GCA and sequence accessions.
    #[serde(rename = "refseq.gcf")]
    RefSeqGcf,
    /// UniProtKB six- or ten-character accession, without isoform/other suffixes.
    #[serde(rename = "uniprot")]
    UniProt,
}

impl Namespace {
    /// Canonical identifiers.org namespace spelling.
    pub const fn prefix(self) -> &'static str {
        match self {
            Self::RefSeqGcf => "refseq.gcf",
            Self::UniProt => "uniprot",
        }
    }
}

impl fmt::Display for Namespace {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.prefix())
    }
}

impl FromStr for Namespace {
    type Err = IdentifierError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value.eq_ignore_ascii_case("refseq.gcf") {
            Ok(Self::RefSeqGcf)
        } else if value.eq_ignore_ascii_case("uniprot") {
            Ok(Self::UniProt)
        } else {
            Err(IdentifierError::UnsupportedNamespace(value.to_owned()))
        }
    }
}

/// Locally validated biological reference; accessions retain case and revision.
///
/// JSON is `{"namespace":"refseq.gcf","accession":"GCF_000001405.40"}`.
/// Deserialization validates the accession, just like [`Self::new`]. Version
/// and rendered URIs are derived, so serialized fields cannot disagree.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(try_from = "IdentifierParts")]
pub struct IdentifierRef {
    namespace: Namespace,
    accession: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct IdentifierParts {
    namespace: Namespace,
    accession: String,
}

impl TryFrom<IdentifierParts> for IdentifierRef {
    type Error = IdentifierError;

    fn try_from(parts: IdentifierParts) -> Result<Self, Self::Error> {
        Self::new(parts.namespace, &parts.accession)
    }
}

impl IdentifierRef {
    /// Validate an accession in an explicitly selected namespace.
    ///
    /// GCF syntax follows the packaged registry: `GCF_`, nine ASCII digits,
    /// optionally a dot and one or more ASCII digits. Numeric suffixes are kept
    /// verbatim; zero and leading zeroes are not proof of an assigned revision.
    /// UniProt syntax follows the provider's six/ten-character accession rules.
    pub fn new(namespace: Namespace, accession: &str) -> Result<Self, IdentifierError> {
        let valid = match namespace {
            Namespace::RefSeqGcf => valid_gcf(accession),
            Namespace::UniProt => valid_uniprot(accession),
        };
        if !valid {
            return Err(IdentifierError::InvalidAccession {
                namespace,
                accession: accession.to_owned(),
            });
        }
        Ok(Self {
            namespace,
            accession: accession.to_owned(),
        })
    }

    /// Parse one complete reference without network I/O or prompt extraction.
    ///
    /// Accepts `namespace:accession`, `namespace://accession`,
    /// `identifiers://namespace:accession`, and HTTP(S) identifiers.org URLs
    /// containing `namespace:accession` or `namespace/accession`. Namespace and
    /// URL scheme/host comparisons are ASCII case-insensitive; accession case
    /// is never changed. Only GCF permits a bare accession. Outer whitespace is
    /// trimmed, but prose, provider prefixes, percent encoding, query strings,
    /// fragments, trailing slashes, and unrecognized namespaces are rejected.
    pub fn parse(value: &str) -> Result<Self, IdentifierError> {
        let value = value.trim();
        if value.is_empty()
            || value.chars().any(|c| c.is_whitespace() || c.is_control())
            || value.contains(['?', '#', '%', '\\'])
        {
            return Err(IdentifierError::InvalidReference);
        }

        if let Some((scheme, rest)) = value.split_once("://") {
            if scheme.eq_ignore_ascii_case("http") || scheme.eq_ignore_ascii_case("https") {
                let (host, path) = rest
                    .split_once('/')
                    .ok_or(IdentifierError::InvalidReference)?;
                if !host.eq_ignore_ascii_case("identifiers.org")
                    && !host.eq_ignore_ascii_case("www.identifiers.org")
                {
                    return Err(IdentifierError::InvalidReference);
                }
                return Self::parse_explicit(path, true);
            }
            if scheme.eq_ignore_ascii_case("identifiers") {
                return Self::parse_explicit(rest, false);
            }
            return Self::new(scheme.parse()?, rest);
        }
        if value.contains(':') {
            return Self::parse_explicit(value, false);
        }
        if value.starts_with("GCF_") {
            return Self::new(Namespace::RefSeqGcf, value);
        }
        Err(IdentifierError::MissingNamespace)
    }

    fn parse_explicit(value: &str, allow_legacy_path: bool) -> Result<Self, IdentifierError> {
        let pair = value.split_once(':').or_else(|| {
            if allow_legacy_path {
                value.split_once('/')
            } else {
                None
            }
        });
        let (namespace, accession) = pair.ok_or(IdentifierError::InvalidReference)?;
        Self::new(namespace.parse()?, accession)
    }

    /// Namespace assigned explicitly or inferred by the narrow bare-GCF rule.
    pub const fn namespace(&self) -> Namespace {
        self.namespace
    }

    /// Full accession, including any assembly revision suffix.
    pub fn accession(&self) -> &str {
        &self.accession
    }

    /// Accession without a GCF revision suffix; unchanged for UniProtKB.
    pub fn base_accession(&self) -> &str {
        self.accession
            .split_once('.')
            .map_or(self.accession.as_str(), |(base, _)| base)
    }

    /// GCF assembly revision text, if explicitly given, preserved exactly.
    ///
    /// `None` does not identify a specific assembly revision. UniProtKB entry
    /// and sequence versions are separate provider metadata, not encoded by the
    /// supported UniProt accession. This method makes no claim of immutability.
    pub fn version(&self) -> Option<&str> {
        self.accession.split_once('.').map(|(_, version)| version)
    }

    /// Canonical compact identifier, with no provider preference.
    pub fn compact_id(&self) -> String {
        format!("{}:{}", self.namespace, self.accession)
    }

    /// Canonical BioV biological data resource URI, never a dataset handle.
    pub fn resource_uri(&self) -> String {
        format!("{}://{}", self.namespace, self.accession)
    }
}

impl FromStr for IdentifierRef {
    type Err = IdentifierError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        Self::parse(value)
    }
}

impl fmt::Display for IdentifierRef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.namespace, self.accession)
    }
}

/// Stable local validation categories. None of these errors reflects network I/O.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IdentifierError {
    /// Empty, decorated, multiple, or otherwise unsupported reference form.
    InvalidReference,
    /// An explicit namespace is required (except for allowlisted bare GCF IDs).
    MissingNamespace,
    /// Namespace is not in this crate's curated subset.
    UnsupportedNamespace(String),
    /// The accession does not match this namespace's supported syntax.
    InvalidAccession {
        namespace: Namespace,
        accession: String,
    },
}

impl fmt::Display for IdentifierError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidReference => f.write_str("expected one supported, undecorated identifier reference"),
            Self::MissingNamespace => f.write_str("an explicit biological namespace is required; only GCF assembly accessions may be bare"),
            Self::UnsupportedNamespace(namespace) => write!(f, "unsupported biological namespace {namespace:?}; supported: refseq.gcf, uniprot"),
            Self::InvalidAccession { namespace, accession } => write!(f, "invalid accession {accession:?} for supported {namespace} syntax"),
        }
    }
}

impl Error for IdentifierError {}

fn valid_gcf(accession: &str) -> bool {
    let (base, version) = accession
        .split_once('.')
        .map_or((accession, None), |(base, version)| (base, Some(version)));
    base.strip_prefix("GCF_")
        .is_some_and(|digits| digits.len() == 9 && digits.bytes().all(|c| c.is_ascii_digit()))
        && version
            .is_none_or(|digits| !digits.is_empty() && digits.bytes().all(|c| c.is_ascii_digit()))
}

fn valid_uniprot(accession: &str) -> bool {
    // Provider rule: [OPQ][0-9][A-Z0-9]{3}[0-9] |
    //               [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
    let bytes = accession.as_bytes();
    if !matches!(bytes.len(), 6 | 10)
        || !bytes
            .iter()
            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit())
        || !bytes[0].is_ascii_uppercase()
        || !bytes[1].is_ascii_digit()
    {
        return false;
    }
    if matches!(bytes[0], b'O' | b'P' | b'Q') {
        return bytes.len() == 6 && bytes[5].is_ascii_digit();
    }
    bytes[2..]
        .chunks_exact(4)
        .all(|chunk| chunk[0].is_ascii_uppercase() && chunk[3].is_ascii_digit())
}
