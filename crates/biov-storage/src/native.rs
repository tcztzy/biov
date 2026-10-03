use crate::{
    tree::{invalid, limit, package, read_bounded, relative},
    EntryKind, InventoryEntry, NativeDeclaration, StorageError, Validation, MAX_REPRESENTATIONS,
    MAX_RESOLVED_FILES,
};
use biov_identifiers::{IdentifierRef, Namespace};
use md5::{Digest, Md5};
use serde::Deserialize;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::Read,
    path::Path,
};

#[derive(Debug, Clone)]
pub(crate) struct Reference {
    pub namespace: &'static str,
    pub accession: String,
    pub base: String,
    pub version: Option<String>,
}
impl Reference {
    pub fn compact(&self) -> String {
        format!("{}:{}", self.namespace, self.accession)
    }
    pub fn matches(&self, canonical: &Self) -> bool {
        self.namespace == canonical.namespace
            && self.base == canonical.base
            && self
                .version
                .as_ref()
                .is_none_or(|v| Some(v) == canonical.version.as_ref())
    }
}

pub(crate) fn parse_ref(value: &str, canonical: bool) -> Result<Reference, StorageError> {
    if value.len() > 128 {
        return Err(invalid("reference exceeds 128 bytes"));
    }
    if let Some((namespace, accession)) = value.split_once(':') {
        if namespace == "pdb" {
            if accession.len() != 4
                || !matches!(accession.as_bytes()[0], b'1'..=b'9')
                || !accession.bytes().all(|v| v.is_ascii_alphanumeric())
            {
                return Err(invalid(
                    "this PDB slice accepts only classic four-character entry accessions",
                ));
            }
            let accession = accession.to_ascii_uppercase();
            return Ok(Reference {
                namespace: "pdb",
                base: accession.clone(),
                accession,
                version: None,
            });
        }
    }
    let id = IdentifierRef::parse(value)
        .map_err(|_| invalid("expected refseq.gcf assembly reference or pdb:entry reference"))?;
    if id.namespace() != Namespace::RefSeqGcf {
        return Err(invalid("unsupported native-store namespace"));
    }
    if canonical && id.version().is_none() {
        return Err(invalid(
            "canonical RefSeq GCF reference must include the assembly version",
        ));
    }
    Ok(Reference {
        namespace: "refseq.gcf",
        accession: id.accession().to_owned(),
        base: id.base_accession().to_owned(),
        version: id.version().map(str::to_owned),
    })
}

pub(crate) fn validate_scope(namespace: &str, scope: &str) -> Result<(), StorageError> {
    if namespace == "refseq.gcf" && scope == "assembly" || namespace == "pdb" && scope == "entry" {
        return Ok(());
    }
    if namespace == "pdb"
        && scope.strip_prefix("assembly:").is_some_and(|id| {
            !id.is_empty()
                && id.len() <= 20
                && !id.starts_with('0')
                && id.bytes().all(|b| b.is_ascii_digit())
        })
    {
        return Ok(());
    }
    Err(invalid(
        "scope must be assembly for RefSeq or entry/assembly:<positive ID> for PDB",
    ))
}

pub(crate) fn representation_name(name: &str) -> Result<(), StorageError> {
    if name.is_empty()
        || name.len() > 64
        || !name
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'_')
    {
        return Err(invalid(
            "representation names must use 1-64 lowercase ASCII letters, digits or underscores",
        ));
    }
    Ok(())
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) struct NativeMetadata {
    pub scope: String,
    pub metadata: Vec<String>,
    pub representations: BTreeMap<String, Vec<String>>,
    pub unavailable: Vec<String>,
    pub validation: Validation,
}

pub(crate) fn validate(
    source: &Path,
    canonical: &Reference,
    declaration: &NativeDeclaration,
    inventory: &[InventoryEntry],
) -> Result<NativeMetadata, StorageError> {
    let files: BTreeMap<&str, &InventoryEntry> = inventory
        .iter()
        .filter(|e| e.kind == EntryKind::File)
        .map(|e| (e.path.as_str(), e))
        .collect();
    match declaration {
        NativeDeclaration::Refseq => {
            if canonical.namespace != "refseq.gcf" {
                return Err(invalid(
                    "RefSeq declaration requires a refseq.gcf reference",
                ));
            }
            refseq(source, canonical, &files)
        }
        NativeDeclaration::Pdb {
            scope,
            representations,
        } => {
            if canonical.namespace != "pdb" {
                return Err(invalid("PDB declaration requires a pdb reference"));
            }
            validate_scope("pdb", scope)?;
            if representations.is_empty() {
                return Err(invalid("PDB declarations need at least one representation"));
            }
            validate_mappings(representations, &files)?;
            let representations = representations.clone();
            Ok(NativeMetadata { scope: scope.clone(), metadata: Vec::new(), representations, unavailable: Vec::new(), validation: Validation {
                method: "caller_declared_pdb_identity_scope_and_representations".into(),
                limits: vec!["Classic PDB accession syntax, relative paths, nonempty mapped files and local byte integrity only".into(),
                    "No mmCIF/PDB parsing, biological accession or assembly verification, provider authenticity, or scientific QC".into(),
                    "Compression and decoding relations are not inferred; native formats remain authoritative".into()],
            } })
        }
    }
}

fn validate_mappings(
    representations: &BTreeMap<String, Vec<String>>,
    files: &BTreeMap<&str, &InventoryEntry>,
) -> Result<(), StorageError> {
    if representations.len() > MAX_REPRESENTATIONS {
        return Err(limit("more than 128 representations"));
    }
    for (name, paths) in representations {
        representation_name(name)?;
        if paths.is_empty() || paths.len() > MAX_RESOLVED_FILES {
            return Err(limit("each representation needs 1-256 files"));
        }
        let mut seen = BTreeSet::new();
        for path in paths {
            relative(path)?;
            if !seen.insert(path) {
                return Err(invalid("duplicate path in representation mapping"));
            }
            if files.get(path.as_str()).is_none_or(|e| e.bytes == 0) {
                return Err(package(
                    "declared representation is missing, empty or not a regular file",
                ));
            }
        }
    }
    Ok(())
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct Catalog {
    api_version: String,
    assemblies: Vec<Assembly>,
}
#[derive(Deserialize)]
struct Assembly {
    accession: Option<String>,
    files: Vec<CatalogFile>,
}
#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct CatalogFile {
    file_path: String,
    file_type: String,
    uncompressed_length_bytes: NativeLength,
}
#[derive(Deserialize)]
#[serde(untagged)]
enum NativeLength {
    String(String),
    Integer(u64),
}
impl NativeLength {
    fn get(&self) -> Result<u64, StorageError> {
        match self {
            Self::Integer(n) => Ok(*n),
            Self::String(n) => n
                .parse()
                .map_err(|_| package("catalog length is not an unsigned byte count")),
        }
    }
}

fn refseq(
    source: &Path,
    canonical: &Reference,
    files: &BTreeMap<&str, &InventoryEntry>,
) -> Result<NativeMetadata, StorageError> {
    const CATALOG: &str = "ncbi_dataset/data/dataset_catalog.json";
    const BASELINE: &[&str] = &[
        "genome_fasta",
        "gff3",
        "protein_fasta",
        "cds_fasta",
        "rna_fasta",
        "assembly_report",
        "sequence_report",
    ];
    if !files.contains_key("README.md")
        || !files.contains_key("md5sum.txt")
        || !files.contains_key(CATALOG)
    {
        return Err(package("RefSeq package requires native README.md, md5sum.txt and ncbi_dataset/data/dataset_catalog.json"));
    }
    let catalog: Catalog = serde_json::from_slice(&read_bounded(&source.join(CATALOG))?)
        .map_err(|_| package("unsupported or malformed native catalog JSON"))?;
    if catalog.api_version != "V2" || catalog.assemblies.is_empty() {
        return Err(package("expected a nonempty V2 NCBI Datasets catalog"));
    }
    let mut required = BTreeSet::from([CATALOG.to_owned()]);
    let mut representations: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut found_accession = false;
    let mut metadata = vec!["README.md".into(), "md5sum.txt".into(), CATALOG.into()];
    for assembly in catalog.assemblies {
        if let Some(accession) = &assembly.accession {
            if accession != &canonical.accession {
                return Err(package("catalog assembly accession does not match the declared canonical versioned GCF"));
            }
            if assembly.files.is_empty() {
                return Err(package(
                    "canonical assembly catalog group has no materialized file declarations",
                ));
            }
            found_accession = true;
        }
        for file in assembly.files {
            if file.file_type.trim().is_empty() {
                return Err(package("catalog fileType must not be empty"));
            }
            if assembly.accession.is_none() && file.file_type != "DATA_REPORT" {
                return Err(package(
                    "accessionless catalog groups may declare only global DATA_REPORT metadata",
                ));
            }
            relative(&file.file_path)?;
            let path = format!("ncbi_dataset/data/{}", file.file_path);
            if !required.insert(path.clone()) {
                return Err(package("duplicate catalog path"));
            }
            let actual = files.get(path.as_str()).ok_or_else(|| package("catalog file is not materialized; dehydrated/incomplete packages are not ready"))?;
            let expected = file.uncompressed_length_bytes.get()?;
            if actual.bytes != expected || expected == 0 {
                return Err(package(
                    "catalog byte length mismatch or empty required native file",
                ));
            }
            let representation = match file.file_type.as_str() {
                "GENOMIC_NUCLEOTIDE_FASTA" => "genome_fasta".to_owned(),
                "RNA_NUCLEOTIDE_FASTA" => "rna_fasta".to_owned(),
                "CDS_NUCLEOTIDE_FASTA" => "cds_fasta".to_owned(),
                "PROTEIN_FASTA" => "protein_fasta".to_owned(),
                "GFF3" => "gff3".to_owned(),
                "DATA_REPORT" => {
                    metadata.push(path.clone());
                    "assembly_report".to_owned()
                }
                "SEQUENCE_REPORT" => {
                    metadata.push(path.clone());
                    "sequence_report".to_owned()
                }
                other => format!("native_{}", other.to_ascii_lowercase()),
            };
            representation_name(&representation)
                .map_err(|_| package("catalog fileType cannot be represented safely"))?;
            representations
                .entry(representation)
                .or_default()
                .push(path);
        }
    }
    if !found_accession || representations.is_empty() {
        return Err(package(
            "catalog does not contain a materialized canonical assembly group",
        ));
    }
    validate_mappings(&representations, files)?;
    let md5_text = String::from_utf8(read_bounded(&source.join("md5sum.txt"))?)
        .map_err(|_| package("native md5sum.txt is not UTF-8"))?;
    let mut checked = BTreeSet::new();
    for line in md5_text.lines() {
        if line.is_empty() {
            continue;
        }
        let raw = line.as_bytes();
        if raw.len() < 35
            || raw[32] != b' '
            || !matches!(raw[33], b' ' | b'*')
            || !raw[..32].iter().all(u8::is_ascii_hexdigit)
        {
            return Err(package("unsupported native MD5 checksum line"));
        }
        let path = &line[34..];
        relative(path)?;
        if !checked.insert(path.to_owned()) {
            return Err(package("duplicate native MD5 checksum path"));
        }
        if !files.contains_key(path) {
            return Err(package("native MD5 target is not materialized"));
        }
        let mut input =
            File::open(source.join(path)).map_err(|e| crate::io("open MD5 target", e))?;
        let mut hash = Md5::new();
        let mut buffer = [0; 64 * 1024];
        loop {
            let n = input
                .read(&mut buffer)
                .map_err(|e| crate::io("verify native MD5", e))?;
            if n == 0 {
                break;
            }
            hash.update(&buffer[..n]);
        }
        if !format!("{:x}", hash.finalize()).eq_ignore_ascii_case(&line[..32]) {
            return Err(package("native MD5 mismatch"));
        }
    }
    if !required.is_subset(&checked) {
        return Err(package(
            "native MD5 inventory must cover the catalog and every catalog file",
        ));
    }
    let unavailable = BASELINE
        .iter()
        .filter(|name| !representations.contains_key(**name))
        .map(|name| (*name).into())
        .collect();
    Ok(NativeMetadata { scope: "assembly".into(), metadata, representations, unavailable, validation: Validation {
        method: "ncbi_datasets_v2_catalog_lengths_and_supplied_md5".into(),
        limits: vec!["The versioned canonical GCF must match every accession-bearing catalog group; global DATA_REPORT groups may omit accession".into(),
            "All catalog paths must be materialized with exact catalog byte lengths and supplied MD5 checksums".into(),
            "Native README, catalog and reports remain authoritative; no sequence/annotation biological QC or source authenticity is established".into(),
            "An absent representation (including RNA) describes local package availability, not biological absence".into()],
    } })
}
