//! Independently anchored syntax examples from provider documentation and the
//! existing Python BioV contract tests. No resolver or provider I/O is needed.

use biov_identifiers::{IdentifierError, IdentifierRef, Namespace};
use serde_json::json;

#[test]
fn python_gcf_reference_forms_share_identity() {
    // tests/test_artifacts.py and tests/test_identifiers_mcp.py use these forms.
    let expected = IdentifierRef::new(Namespace::RefSeqGcf, "GCF_000001030.2").unwrap();
    for value in [
        "GCF_000001030.2",
        "refseq.gcf:GCF_000001030.2",
        "refseq.gcf://GCF_000001030.2",
        "identifiers://refseq.gcf:GCF_000001030.2",
        "https://identifiers.org/refseq.gcf:GCF_000001030.2",
        "http://www.identifiers.org/refseq.gcf/GCF_000001030.2",
        " HTTPS://IDENTIFIERS.ORG/REFSEQ.GCF:GCF_000001030.2 \n",
    ] {
        assert_eq!(IdentifierRef::parse(value).unwrap(), expected, "{value}");
    }
    assert_eq!(expected.compact_id(), "refseq.gcf:GCF_000001030.2");
    assert_eq!(expected.resource_uri(), "refseq.gcf://GCF_000001030.2");
    assert_eq!(expected.base_accession(), "GCF_000001030");
    assert_eq!(expected.version(), Some("2"));
}

#[test]
fn python_uniprot_requires_an_explicit_namespace() {
    for value in [
        "uniprot:P42212",
        "UNIPROT://P42212",
        "identifiers://uniprot:P42212",
        "https://identifiers.org/uniprot/P42212",
    ] {
        let parsed = IdentifierRef::parse(value).unwrap();
        assert_eq!(parsed.namespace(), Namespace::UniProt);
        assert_eq!(parsed.accession(), "P42212");
        assert_eq!(parsed.base_accession(), "P42212");
        assert_eq!(parsed.version(), None);
    }
    assert_eq!(
        IdentifierRef::parse("P12345"),
        Err(IdentifierError::MissingNamespace)
    );
}

#[test]
fn official_uniprot_examples_validate() {
    // https://www.uniprot.org/help/accession_numbers lists these exact examples.
    for accession in ["A2BC19", "P12345", "A0A023GPI8"] {
        assert!(IdentifierRef::new(Namespace::UniProt, accession).is_ok());
    }
    // Packaged identifiers.org registry sample: src/biov/assets/...registry.json.
    assert!(IdentifierRef::new(Namespace::UniProt, "P0DP23").is_ok());
}

#[test]
fn uniprot_full_pattern_rejects_wrong_positions_and_non_ascii() {
    // Unlike the packaged upstream regex, the curated provider rule contains
    // no comma or space in accession character classes.
    for accession in [
        "P1234",
        "P123456",
        "P1234A",
        "P1,,A5",
        "P1  A5",
        "p12345",
        "P１２３４５",
        "12BC19",
        "ABBC19",
        "A21C19",
        "A2BC1Z",
        "A0A0239PI8",
        "A0A023GPIX",
        "P0A023GPI8",
        "Q0A023GPI8",
        "O0A023GPI8",
    ] {
        assert!(
            IdentifierRef::new(Namespace::UniProt, accession).is_err(),
            "unexpected valid accession {accession}"
        );
    }
}

#[test]
fn uniprot_isoforms_revisions_and_fragments_are_explicitly_out_of_scope() {
    for value in [
        "uniprot:P12345-2",
        "uniprot:P12345.1",
        "uniprot:P12345#PRO_000001",
        "uniprot:A0A023GPI8-1",
        "uniprot:EGFR_HUMAN",
    ] {
        assert!(IdentifierRef::parse(value).is_err(), "{value}");
    }
}

#[test]
fn official_ncbi_refseq_assembly_example_preserves_revision() {
    // NCBI genome-processing/version-status documents GCF_009914755.1 and its
    // separately versioned GenBank partner GCA_009914755.4.
    let assembly = IdentifierRef::parse("GCF_009914755.1").unwrap();
    assert_eq!(assembly.version(), Some("1"));
    assert!(IdentifierRef::new(Namespace::RefSeqGcf, "GCA_009914755.4").is_err());
}

#[test]
fn gcf_version_absence_and_differences_are_not_collapsed() {
    let bare = IdentifierRef::parse("GCF_000001405").unwrap();
    let first = IdentifierRef::parse("GCF_000001405.1").unwrap();
    let second = IdentifierRef::parse("GCF_000001405.2").unwrap();
    assert_eq!(bare.version(), None);
    assert_ne!(bare, first);
    assert_ne!(first, second);
    assert_eq!(first.base_accession(), second.base_accession());
}

#[test]
fn gcf_keeps_numeric_revision_lexemes_without_asserting_they_exist() {
    // The packaged GCF registry pattern is ^GCF_[0-9]{9}(\.[0-9]+)?$.
    // Do not silently normalize leading zeroes, reject by integer overflow,
    // or infer that a syntactically valid revision has been assigned by NCBI.
    for version in ["0", "01", "1234567890123456789012345678901234567890"] {
        let accession = format!("GCF_000001405.{version}");
        let reference = IdentifierRef::parse(&accession).unwrap();
        assert_eq!(reference.version(), Some(version));
        assert_eq!(reference.accession(), accession);
    }
}

#[test]
fn invalid_gcf_envelopes_reject_before_io() {
    for value in [
        "GCF_00001030", // Existing Python rejection: only eight digits.
        "GCF_0000010300",
        "GCF_000001030.",
        "GCF_000001030.1.2",
        "GCF_000001030.-2",
        "GCF_000001030.+2",
        "GCF_000001030.２",
        "GCF_00000103A.2",
        "gcf_000001030.2",
        "NC_000001.11",
    ] {
        assert!(
            IdentifierRef::new(Namespace::RefSeqGcf, value).is_err(),
            "{value}"
        );
    }
}

#[test]
fn whole_reference_contract_rejects_prose_and_multiple_values() {
    for value in [
        "",
        "   ",
        "What is the GC content of GCF_000001030.2?",
        "GCF_000001030.2 GCF_000001030.2",
        "refseq.gcf:GCF_000001030.2,",
        "uniprot:P12345\nuniprot:P42212",
        "uniprot:P12\0A5",
    ] {
        assert!(IdentifierRef::parse(value).is_err(), "{value:?}");
    }
}

#[test]
fn urls_are_exact_and_do_not_discard_scientific_qualifiers() {
    for value in [
        "https://identifiers.org/uniprot:P12345?version=2",
        "https://identifiers.org/uniprot:P12345#PRO_000001",
        "https://identifiers.org/uniprot:P12345/",
        "https://identifiers.org/uniprot:P%31%32%33%34%35",
        "https://identifiers.org.evil.test/uniprot:P12345",
        "https://identifiers.org@evil.test/uniprot:P12345",
        "https://example.org/uniprot:P12345",
        "https://identifiers.org",
        "https://identifiers.org/",
        "https://identifiers.org//uniprot:P12345",
        "https://identifiers.org/ebi/uniprot:P12345",
        "ebi/uniprot:P12345",
        "identifiers://uniprot/P12345",
    ] {
        assert!(IdentifierRef::parse(value).is_err(), "{value}");
    }
}

#[test]
fn dataset_handles_and_paths_are_not_biological_identifiers() {
    for value in [
        "dataset://abc123",
        "biov://datasets/abc123",
        "dataset_abc123",
        "abc123",
        "/tmp/GCF_000001030.2",
        "file:///tmp/table.parquet",
        "doi:10.1038/example",
    ] {
        assert!(IdentifierRef::parse(value).is_err(), "{value}");
    }
}

#[test]
fn serde_roundtrip_contains_only_validated_identity_coordinates() {
    let reference = IdentifierRef::parse("GCF_000001405.40").unwrap();
    let encoded = serde_json::to_value(&reference).unwrap();
    assert_eq!(
        encoded,
        json!({
            "namespace": "refseq.gcf", "accession": "GCF_000001405.40"
        })
    );
    assert_eq!(
        serde_json::from_value::<IdentifierRef>(encoded).unwrap(),
        reference
    );
}

#[test]
fn serde_does_not_bypass_namespace_or_accession_validation() {
    for encoded in [
        json!({"namespace": "uniprot", "accession": "GCF_000001405.40"}),
        json!({"namespace": "refseq.gcf", "accession": "P12345"}),
        json!({"namespace": "uniprot", "accession": "P12345-2"}),
        json!({"namespace": "dataset", "accession": "abc123"}),
        json!({"namespace": "uniprot", "accession": "P12345", "version": 2}),
        json!({"namespace": "uniprot", "accession": "P12345", "resource_uri": "dataset://a"}),
    ] {
        assert!(serde_json::from_value::<IdentifierRef>(encoded).is_err());
    }
}

#[test]
fn display_and_from_str_roundtrip_to_canonical_compact_id() {
    let parsed: IdentifierRef = "UNIPROT://P12345".parse().unwrap();
    assert_eq!(parsed.to_string(), "uniprot:P12345");
    assert_eq!(parsed.to_string().parse::<IdentifierRef>().unwrap(), parsed);
    assert_eq!(
        "REFSEQ.GCF".parse::<Namespace>().unwrap().prefix(),
        "refseq.gcf"
    );
}
