use crate::*;
use md5::{Digest, Md5};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
    sync::{Arc, Barrier},
};
use tempfile::TempDir;

struct Fixture {
    _temp: TempDir,
    source_root: PathBuf,
    store_root: PathBuf,
}
impl Fixture {
    fn new() -> Self {
        let temp = tempfile::tempdir().unwrap();
        let source_root = temp.path().join("imports");
        let store_root = temp.path().join("store");
        fs::create_dir(&source_root).unwrap();
        fs::create_dir(&store_root).unwrap();
        Self {
            _temp: temp,
            source_root,
            store_root,
        }
    }
    fn store(&self) -> NativeStore {
        NativeStore::new(&self.store_root).unwrap()
    }
    fn refseq(&self) -> RegisterRequest {
        self.refseq_version("GCF_000005845.2", "package")
    }
    fn refseq_version(&self, accession: &str, folder: &str) -> RegisterRequest {
        let root = self.source_root.join(folder);
        let data = root.join("ncbi_dataset/data");
        let assembly = data.join(accession);
        fs::create_dir_all(&assembly).unwrap();
        fs::create_dir(root.join("empty")).unwrap();
        fs::write(
            root.join("README.md"),
            "Native NCBI-style fixture: FASTA and matching GFF3.\n",
        )
        .unwrap();
        fs::write(assembly.join("genome.fna"), ">chr\nACGT\n").unwrap();
        fs::write(
            assembly.join("genomic.gff"),
            "##gff-version 3\nchr\ttest\tgene\t1\t4\t.\t+\t.\tID=g1\n",
        )
        .unwrap();
        fs::write(
            data.join("assembly_data_report.jsonl"),
            format!("{{\"accession\":\"{accession}\"}}\n"),
        )
        .unwrap();
        write_catalog_and_md5(&root, accession);
        RegisterRequest {
            source_path: folder.into(),
            requested_ref: "refseq.gcf:GCF_000005845".into(),
            canonical_ref: format!("refseq.gcf:{accession}"),
            declaration: NativeDeclaration::Refseq,
        }
    }
    fn pdb(&self) -> RegisterRequest {
        let root = self.source_root.join("pdb");
        fs::create_dir(&root).unwrap();
        // Deliberately not a verified coordinate record: this adapter's limited
        // contract is caller declaration and saved-byte consistency only.
        fs::write(
            root.join("declared.cif"),
            "caller-declared structure bytes\n",
        )
        .unwrap();
        fs::write(root.join("other.cif"), "second caller-declared structure\n").unwrap();
        RegisterRequest {
            source_path: "pdb".into(),
            requested_ref: "pdb:1crn".into(),
            canonical_ref: "pdb:1CRN".into(),
            declaration: NativeDeclaration::Pdb {
                scope: "entry".into(),
                representations: BTreeMap::from([(
                    "structure_cif".into(),
                    vec!["other.cif".into(), "declared.cif".into()],
                )]),
            },
        }
    }
    fn resolve(&self, representation: &str) -> Resolution {
        self.store()
            .resolve(resolve_refseq(representation))
            .unwrap()
    }
}
fn write_catalog_and_md5(root: &Path, accession: &str) {
    let data = root.join("ncbi_dataset/data");
    let genome = format!("{accession}/genome.fna");
    let annotation = format!("{accession}/genomic.gff");
    let catalog = serde_json::json!({"apiVersion":"V2","assemblies":[
        {"files":[{"filePath":"assembly_data_report.jsonl","fileType":"DATA_REPORT","uncompressedLengthBytes":fs::metadata(data.join("assembly_data_report.jsonl")).unwrap().len().to_string()}]},
        {"accession":accession,"files":[
            {"filePath":genome,"fileType":"GENOMIC_NUCLEOTIDE_FASTA","uncompressedLengthBytes":fs::metadata(data.join(&genome)).unwrap().len().to_string()},
            {"filePath":annotation,"fileType":"GFF3","uncompressedLengthBytes":fs::metadata(data.join(&annotation)).unwrap().len().to_string()}]}]});
    fs::write(
        data.join("dataset_catalog.json"),
        serde_json::to_vec_pretty(&catalog).unwrap(),
    )
    .unwrap();
    let mut md5 = String::new();
    for rel in [
        "dataset_catalog.json".to_owned(),
        "assembly_data_report.jsonl".to_owned(),
        genome,
        annotation,
    ] {
        let bytes = fs::read(data.join(&rel)).unwrap();
        md5.push_str(&format!(
            "{:x}  ncbi_dataset/data/{rel}\n",
            Md5::digest(bytes)
        ));
    }
    fs::write(root.join("md5sum.txt"), md5).unwrap();
}
fn resolve_refseq(representation: &str) -> ResolveRequest {
    ResolveRequest {
        reference: "refseq.gcf:GCF_000005845.2".into(),
        representation: representation.into(),
        snapshot_id: None,
        scope: None,
    }
}
fn receipt(f: &Fixture, registration: &Registration) -> Receipt {
    serde_json::from_slice(
        &fs::read(f.store_root.join(&registration.snapshot.receipt_path)).unwrap(),
    )
    .unwrap()
}

#[test]
fn construction_and_missing_resolution_are_read_only() {
    let f = Fixture::new();
    assert!(matches!(f.resolve("genome_fasta"), Resolution::Miss));
    assert_eq!(fs::read_dir(&f.store_root).unwrap().count(), 0);
}
#[test]
fn refseq_copy_is_complete_portable_and_native_metadata_is_authoritative() {
    let f = Fixture::new();
    let request = f.refseq();
    let before = crate::tree::inventory(&f.source_root.join("package"), None).unwrap();
    let result = f.store().register(&f.source_root, request).unwrap();
    assert!(!result.reused);
    let r = receipt(&f, &result);
    assert_eq!(r.inventory, before);
    assert_eq!(
        r.inventory,
        crate::tree::inventory(&f.source_root.join("package"), None).unwrap()
    );
    assert!(r.acquired_at.is_none() && r.source_url.is_none() && r.acquisition_tool.is_none());
    assert!(r.unavailable_representations.contains(&"rna_fasta".into()));
    assert!(r
        .inventory
        .iter()
        .any(|e| e.path == "empty" && e.kind == EntryKind::Directory));
    let Resolution::Ready { snapshot, paths } = f.resolve("genome_fasta") else {
        panic!("expected ready")
    };
    assert_eq!(snapshot, result.snapshot);
    assert_eq!(paths.len(), 1);
    assert_eq!(
        fs::read(&paths[0].execution_host_path).unwrap(),
        b">chr\nACGT\n"
    );
    assert_eq!(
        f.store_root.join(&paths[0].relative_path),
        Path::new(&paths[0].execution_host_path)
    );
    let wrapper =
        fs::read_to_string(f.store_root.join(&snapshot.snapshot_path).join("README.md")).unwrap();
    assert!(
        wrapper.contains("genome_fasta")
            && wrapper.contains("genomic.gff")
            && wrapper.contains("sha256sum -c")
    );
}
#[test]
fn missing_representation_is_unavailable_not_biological_absence() {
    let f = Fixture::new();
    f.store().register(&f.source_root, f.refseq()).unwrap();
    let Resolution::Unavailable {
        available_representations,
        ..
    } = f.resolve("rna_fasta")
    else {
        panic!("expected unavailable")
    };
    assert!(available_representations.contains(&"genome_fasta".into()));
}
#[test]
fn repeated_identical_content_reuses_without_rewriting_first_receipt() {
    let f = Fixture::new();
    let mut request = f.refseq();
    let first = f.store().register(&f.source_root, request.clone()).unwrap();
    let before = fs::read(f.store_root.join(&first.snapshot.receipt_path)).unwrap();
    request.requested_ref = request.canonical_ref.clone();
    let again = f.store().register(&f.source_root, request).unwrap();
    assert!(again.reused);
    assert_eq!(first.snapshot, again.snapshot);
    assert_eq!(again.requested_ref, "refseq.gcf:GCF_000005845.2");
    assert_eq!(
        fs::read(f.store_root.join(&first.snapshot.receipt_path)).unwrap(),
        before
    );
}
#[test]
fn changed_annotation_same_accession_creates_ambiguity_and_explicit_pin_resolves() {
    let f = Fixture::new();
    let request = f.refseq();
    let first = f.store().register(&f.source_root, request.clone()).unwrap();
    fs::write(
        f.source_root
            .join("package/ncbi_dataset/data/GCF_000005845.2/genomic.gff"),
        "##gff-version 3\nchr\ttest\tgene\t1\t3\t.\t+\t.\tID=g2\n",
    )
    .unwrap();
    write_catalog_and_md5(&f.source_root.join("package"), "GCF_000005845.2");
    let second = f.store().register(&f.source_root, request).unwrap();
    assert_ne!(first.snapshot.snapshot_id, second.snapshot.snapshot_id);
    let Resolution::Ambiguous { candidates } = f.resolve("genome_fasta") else {
        panic!("expected ambiguity")
    };
    assert_eq!(candidates.len(), 2);
    let mut resolve = resolve_refseq("genome_fasta");
    resolve.snapshot_id = Some(first.snapshot.snapshot_id.clone());
    let Resolution::Ready { snapshot, .. } = f.store().resolve(resolve).unwrap() else {
        panic!("pin failed")
    };
    assert_eq!(snapshot, first.snapshot);
}
#[test]
fn biological_versions_are_separate_and_unversioned_is_not_latest() {
    let f = Fixture::new();
    let a = f.refseq_version("GCF_000005845.2", "v2");
    let b = f.refseq_version("GCF_000005845.3", "v3");
    f.store().register(&f.source_root, a).unwrap();
    f.store().register(&f.source_root, b).unwrap();
    let mut request = resolve_refseq("genome_fasta");
    request.reference = "refseq.gcf:GCF_000005845".into();
    assert!(matches!(
        f.store().resolve(request).unwrap(),
        Resolution::Ambiguous { .. }
    ));
    assert!(matches!(
        f.resolve("genome_fasta"),
        Resolution::Ready { .. }
    ));
}
#[test]
fn moved_store_and_removed_originals_resolve() {
    let f = Fixture::new();
    f.store().register(&f.source_root, f.refseq()).unwrap();
    let moved = f._temp.path().join("moved");
    fs::rename(&f.store_root, &moved).unwrap();
    fs::remove_dir_all(&f.source_root).unwrap();
    let Resolution::Ready { paths, .. } = NativeStore::new(&moved)
        .unwrap()
        .resolve(resolve_refseq("genome_fasta"))
        .unwrap()
    else {
        panic!("moved store failed")
    };
    assert_eq!(
        fs::read_to_string(&paths[0].execution_host_path).unwrap(),
        ">chr\nACGT\n"
    );
}
#[test]
fn corruption_is_never_replaced_by_reregistering_good_originals() {
    let f = Fixture::new();
    let request = f.refseq();
    let result = f.store().register(&f.source_root, request.clone()).unwrap();
    let file = f
        .store_root
        .join(&result.snapshot.snapshot_path)
        .join("source/ncbi_dataset/data/GCF_000005845.2/genome.fna");
    fs::write(&file, b"damaged").unwrap();
    assert!(matches!(
        f.resolve("genome_fasta"),
        Resolution::Corrupt { .. }
    ));
    assert!(f.store().register(&f.source_root, request).is_err());
    assert_eq!(fs::read(file).unwrap(), b"damaged");
    assert_eq!(
        fs::read(
            f.source_root
                .join("package/ncbi_dataset/data/GCF_000005845.2/genome.fna")
        )
        .unwrap(),
        b">chr\nACGT\n"
    );
}
#[test]
fn wrapper_or_receipt_tampering_and_extra_files_are_corrupt() {
    for target in [
        "README.md",
        "checksums.sha256",
        "acquisition.json",
        "source/unregistered.txt",
    ] {
        let f = Fixture::new();
        let result = f.store().register(&f.source_root, f.refseq()).unwrap();
        fs::write(
            f.store_root
                .join(result.snapshot.snapshot_path)
                .join(target),
            b"changed",
        )
        .unwrap();
        assert!(
            matches!(f.resolve("genome_fasta"), Resolution::Corrupt { .. }),
            "{target}"
        );
    }
}
#[test]
fn pdb_is_explicitly_declared_and_preserves_multiple_paths() {
    let f = Fixture::new();
    let request = f.pdb();
    let result = f.store().register(&f.source_root, request).unwrap();
    assert_eq!(result.snapshot.canonical_ref, "pdb:1CRN");
    let r = receipt(&f, &result);
    assert!(r.validation.method.starts_with("caller_declared"));
    let resolve = ResolveRequest {
        reference: "pdb:1crn".into(),
        representation: "structure_cif".into(),
        snapshot_id: None,
        scope: None,
    };
    let Resolution::Ready { paths, .. } = f.store().resolve(resolve).unwrap() else {
        panic!("expected ready")
    };
    assert_eq!(paths.len(), 2);
}
#[test]
fn same_bytes_with_different_pdb_scope_or_mapping_conflict_without_mutation() {
    for scope_change in [true, false] {
        let f = Fixture::new();
        let request = f.pdb();
        let first = f.store().register(&f.source_root, request.clone()).unwrap();
        let before = fs::read(f.store_root.join(&first.snapshot.receipt_path)).unwrap();
        let mut changed = request;
        let NativeDeclaration::Pdb {
            scope,
            representations,
        } = &mut changed.declaration
        else {
            unreachable!()
        };
        if scope_change {
            *scope = "assembly:1".into();
        } else {
            representations.remove("structure_cif");
            representations.insert("other".into(), vec!["declared.cif".into()]);
        }
        assert!(matches!(
            f.store().register(&f.source_root, changed),
            Err(StorageError::DeclarationConflict)
        ));
        assert_eq!(
            fs::read(f.store_root.join(first.snapshot.receipt_path)).unwrap(),
            before
        );
    }
}
#[test]
fn declared_assembly_never_resolves_as_entry_by_default() {
    let f = Fixture::new();
    let mut request = f.pdb();
    if let NativeDeclaration::Pdb { scope, .. } = &mut request.declaration {
        *scope = "assembly:1".into();
    }
    f.store().register(&f.source_root, request).unwrap();
    let mut resolve = ResolveRequest {
        reference: "pdb:1CRN".into(),
        representation: "structure_cif".into(),
        snapshot_id: None,
        scope: None,
    };
    assert!(matches!(
        f.store().resolve(resolve.clone()).unwrap(),
        Resolution::Miss
    ));
    resolve.scope = Some("assembly:1".into());
    assert!(matches!(
        f.store().resolve(resolve).unwrap(),
        Resolution::Ready { .. }
    ));
}
#[test]
fn dehydrated_wrong_catalog_md5_and_missing_coverage_never_publish() {
    for mode in [
        "missing",
        "length",
        "md5",
        "coverage",
        "canonical",
        "traversal",
    ] {
        let f = Fixture::new();
        let request = f.refseq();
        let root = f.source_root.join("package");
        let genome = root.join("ncbi_dataset/data/GCF_000005845.2/genome.fna");
        match mode {
            "missing" => fs::remove_file(&genome).unwrap(),
            "length" => fs::write(&genome, b">chr\nACGTT\n").unwrap(),
            "md5" => fs::write(&genome, b">chr\nTGCA\n").unwrap(),
            "coverage" => fs::write(root.join("md5sum.txt"), "").unwrap(),
            "canonical" => {
                let path = root.join("ncbi_dataset/data/dataset_catalog.json");
                let text = fs::read_to_string(&path)
                    .unwrap()
                    .replace("GCF_000005845.2", "GCF_000005845.3");
                fs::write(path, text).unwrap();
            }
            "traversal" => {
                let path = root.join("ncbi_dataset/data/dataset_catalog.json");
                let text = fs::read_to_string(&path)
                    .unwrap()
                    .replace("GCF_000005845.2/genome.fna", "../../outside");
                fs::write(path, text).unwrap();
            }
            _ => unreachable!(),
        }
        assert!(
            f.store().register(&f.source_root, request).is_err(),
            "{mode}"
        );
        assert!(matches!(f.resolve("genome_fasta"), Resolution::Miss));
        assert_eq!(
            fs::read_dir(f.store_root.join(".staging")).unwrap().count(),
            0
        );
    }
}
#[test]
fn all_source_and_mapping_traversal_and_nonportable_paths_are_rejected() {
    for path in [
        "../package",
        "/tmp/package",
        "package/../package",
        "package\\child",
        "package\nchild",
        "package//child",
        "package/./child",
    ] {
        let f = Fixture::new();
        let mut request = f.refseq();
        request.source_path = path.into();
        assert!(
            matches!(
                f.store().register(&f.source_root, request),
                Err(StorageError::InvalidInput(_))
            ),
            "{path}"
        );
    }
    let f = Fixture::new();
    let mut request = f.pdb();
    if let NativeDeclaration::Pdb {
        representations, ..
    } = &mut request.declaration
    {
        representations.insert("escape".into(), vec!["../other".into()]);
    }
    assert!(f.store().register(&f.source_root, request).is_err());
}
#[test]
fn overlapping_roots_are_rejected() {
    let f = Fixture::new();
    let request = f.refseq();
    assert!(matches!(
        f.store().register(f._temp.path(), request),
        Err(StorageError::InvalidInput(_))
    ));
}
#[cfg(unix)]
#[test]
fn symlinks_and_special_files_are_rejected() {
    use std::os::unix::fs::symlink;
    for mode in [
        "file",
        "directory",
        "dangling",
        "fifo",
        "source_root",
        "store_root",
    ] {
        let f = Fixture::new();
        let request = f.refseq();
        let package = f.source_root.join("package");
        match mode {
            "file" => symlink(package.join("README.md"), package.join("link")).unwrap(),
            "directory" => symlink(package.join("ncbi_dataset"), package.join("link")).unwrap(),
            "dangling" => symlink(package.join("nonexistent"), package.join("link")).unwrap(),
            "fifo" => {
                rustix::fs::mkfifoat(
                    rustix::fs::CWD,
                    package.join("fifo"),
                    rustix::fs::Mode::RUSR | rustix::fs::Mode::WUSR,
                )
                .unwrap();
            }
            "source_root" => {
                let alias = f._temp.path().join("source_alias");
                symlink(&f.source_root, &alias).unwrap();
                assert!(f.store().register(alias, request).is_err());
                continue;
            }
            "store_root" => {
                let alias = f._temp.path().join("store_alias");
                symlink(&f.store_root, &alias).unwrap();
                assert!(NativeStore::new(alias).is_err());
                continue;
            }
            _ => unreachable!(),
        }
        assert!(
            f.store().register(&f.source_root, request).is_err(),
            "{mode}"
        );
    }
}
#[test]
fn partial_staging_is_never_ready_and_existing_incomplete_destination_is_not_replaced() {
    let f = Fixture::new();
    let request = f.refseq();
    fs::create_dir_all(f.store_root.join(".staging/crashed/source")).unwrap();
    fs::write(
        f.store_root.join(".staging/crashed/acquisition.json"),
        "partial",
    )
    .unwrap();
    assert!(matches!(f.resolve("genome_fasta"), Resolution::Miss));
    let result = f.store().register(&f.source_root, request.clone()).unwrap();
    let snapshot = f.store_root.join(&result.snapshot.snapshot_path);
    fs::remove_dir_all(&snapshot).unwrap();
    fs::create_dir(&snapshot).unwrap();
    assert!(f.store().register(&f.source_root, request).is_err());
    assert_eq!(fs::read_dir(snapshot).unwrap().count(), 0);
    assert!(matches!(
        f.resolve("genome_fasta"),
        Resolution::Corrupt { .. }
    ));
    assert!(f.store_root.join(".staging/crashed").is_dir());
}
#[test]
fn concurrent_identical_registration_has_one_verified_winner() {
    let f = Fixture::new();
    let request = f.refseq();
    let barrier = Arc::new(Barrier::new(6));
    let workers: Vec<_> = (0..6)
        .map(|_| {
            let source = f.source_root.clone();
            let store = f.store_root.clone();
            let request = request.clone();
            let barrier = barrier.clone();
            std::thread::spawn(move || {
                barrier.wait();
                NativeStore::new(store)
                    .unwrap()
                    .register(source, request)
                    .unwrap()
            })
        })
        .collect();
    let results: Vec<_> = workers.into_iter().map(|w| w.join().unwrap()).collect();
    assert_eq!(results.iter().filter(|r| !r.reused).count(), 1);
    assert!(results.iter().all(|r| r.snapshot == results[0].snapshot));
    assert_eq!(
        fs::read_dir(
            f.store_root
                .join("artifacts/refseq.gcf/GCF_000005845.2/snapshots")
        )
        .unwrap()
        .count(),
        1
    );
    assert_eq!(
        fs::read_dir(f.store_root.join(".staging")).unwrap().count(),
        0
    );
}
#[test]
fn empty_directories_participate_in_source_identity() {
    let f = Fixture::new();
    let request = f.refseq();
    let first = f.store().register(&f.source_root, request.clone()).unwrap();
    fs::create_dir(f.source_root.join("package/another_empty")).unwrap();
    let second = f.store().register(&f.source_root, request).unwrap();
    assert_ne!(first.snapshot.snapshot_id, second.snapshot.snapshot_id);
}
#[test]
fn streamed_native_files_are_not_subject_to_table_64_mib_limit() {
    let f = Fixture::new();
    let request = f.pdb();
    let file = fs::OpenOptions::new()
        .write(true)
        .open(f.source_root.join("pdb/declared.cif"))
        .unwrap();
    file.set_len(65 * 1024 * 1024).unwrap();
    drop(file);
    let registration = f.store().register(&f.source_root, request).unwrap();
    assert!(registration.snapshot.total_bytes > 64 * 1024 * 1024);
    let resolution = f
        .store()
        .resolve(ResolveRequest {
            reference: "pdb:1CRN".into(),
            representation: "structure_cif".into(),
            snapshot_id: None,
            scope: None,
        })
        .unwrap();
    assert!(matches!(resolution, Resolution::Ready { .. }));
}
#[test]
fn metadata_count_and_request_limits_are_bounded() {
    let f = Fixture::new();
    let mut request = f.pdb();
    if let NativeDeclaration::Pdb {
        representations, ..
    } = &mut request.declaration
    {
        for n in 0..MAX_REPRESENTATIONS {
            representations.insert(format!("r{n}"), vec!["declared.cif".into()]);
        }
    }
    assert!(matches!(
        f.store().register(&f.source_root, request),
        Err(StorageError::Limit(_))
    ));
    let mut request = resolve_refseq("genome_fasta");
    request.snapshot_id = Some("sha256-../../x".into());
    assert!(f.store().resolve(request).is_err());
    let path = f.source_root.join("oversize.json");
    let file = fs::File::create(&path).unwrap();
    file.set_len(MAX_METADATA_BYTES as u64 + 1).unwrap();
    assert!(matches!(
        crate::tree::read_bounded(&path),
        Err(StorageError::Limit(_))
    ));
}
#[test]
fn source_inventory_comparison_detects_addition_removal_and_same_size_edit() {
    let f = Fixture::new();
    f.pdb();
    let source = f.source_root.join("pdb");
    let copy = f.store_root.join("temporary");
    fs::create_dir(&copy).unwrap();
    let before = crate::tree::inventory(&source, Some(&copy)).unwrap();
    fs::write(
        source.join("declared.cif"),
        "Caller-declared structure bytes\n",
    )
    .unwrap();
    assert_ne!(before, crate::tree::inventory(&source, None).unwrap());
    fs::write(
        source.join("declared.cif"),
        "caller-declared structure bytes\n",
    )
    .unwrap();
    fs::write(source.join("new"), "x").unwrap();
    assert_ne!(before, crate::tree::inventory(&source, None).unwrap());
    fs::remove_file(source.join("new")).unwrap();
    fs::remove_file(source.join("other.cif")).unwrap();
    assert_ne!(before, crate::tree::inventory(&source, None).unwrap());
}

#[test]
fn inventory_encoding_matches_independent_python_hashlib_golden() {
    let entries = vec![
        InventoryEntry {
            path: "empty".into(),
            kind: EntryKind::Directory,
            bytes: 0,
            sha256: None,
        },
        InventoryEntry {
            path: "x.txt".into(),
            kind: EntryKind::File,
            bytes: 3,
            sha256: Some("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad".into()),
        },
    ];
    // Independently calculated with Python hashlib + struct.pack('>I' / '>Q').
    assert_eq!(
        crate::tree::content_digest(&entries).unwrap(),
        "b12e290941cc08d8e2aa5ebcae6c9787cc0c5f4dcd77708dfd9723f3554ec868"
    );
}
#[test]
fn representation_order_is_preserved_and_reordered_declaration_conflicts() {
    let f = Fixture::new();
    let mut request = f.pdb();
    let first = f.store().register(&f.source_root, request.clone()).unwrap();
    assert_eq!(
        receipt(&f, &first).representations["structure_cif"],
        ["other.cif", "declared.cif"]
    );
    if let NativeDeclaration::Pdb {
        representations, ..
    } = &mut request.declaration
    {
        representations.get_mut("structure_cif").unwrap().reverse();
    }
    assert!(matches!(
        f.store().register(&f.source_root, request),
        Err(StorageError::DeclarationConflict)
    ));
}
#[test]
fn root_readme_is_checked_before_publish_and_before_reuse() {
    let f = Fixture::new();
    let request = f.refseq();
    fs::create_dir(f.store_root.join("README.md")).unwrap();
    assert!(f.store().register(&f.source_root, request.clone()).is_err());
    assert!(matches!(f.resolve("genome_fasta"), Resolution::Miss));
    fs::remove_dir(f.store_root.join("README.md")).unwrap();
    f.store().register(&f.source_root, request.clone()).unwrap();
    fs::remove_file(f.store_root.join("README.md")).unwrap();
    fs::create_dir(f.store_root.join("README.md")).unwrap();
    assert!(f.store().register(&f.source_root, request).is_err());
}
#[test]
fn canonical_and_pinned_lookups_do_not_scan_unrelated_directory_entries() {
    let f = Fixture::new();
    let result = f.store().register(&f.source_root, f.refseq()).unwrap();
    let namespace = f.store_root.join("artifacts/refseq.gcf");
    for n in 0..=MAX_SCAN_CANDIDATES {
        fs::create_dir(namespace.join(format!("unrelated-{n}"))).unwrap();
    }
    assert!(matches!(
        f.resolve("genome_fasta"),
        Resolution::Ready { .. }
    ));
    let mut request = resolve_refseq("genome_fasta");
    request.reference = "refseq.gcf:GCF_000005845".into();
    assert!(matches!(
        f.store().resolve(request),
        Err(StorageError::Limit(_))
    ));
    let snapshots = f
        .store_root
        .join("artifacts/refseq.gcf/GCF_000005845.2/snapshots");
    for n in 0..=MAX_SCAN_CANDIDATES {
        fs::create_dir(snapshots.join(format!("ignored-{n}"))).unwrap();
    }
    assert!(matches!(
        f.store().resolve(resolve_refseq("genome_fasta")),
        Err(StorageError::Limit(_))
    ));
    let mut pinned = resolve_refseq("genome_fasta");
    pinned.snapshot_id = Some(result.snapshot.snapshot_id);
    assert!(matches!(
        f.store().resolve(pinned).unwrap(),
        Resolution::Ready { .. }
    ));
}
#[cfg(unix)]
#[test]
fn dangling_store_symlinks_never_create_external_files() {
    use std::os::unix::fs::symlink;
    let f = Fixture::new();
    let request = f.refseq();
    fs::create_dir(f.store_root.join(".locks")).unwrap();
    let outside = f.source_root.join("must_not_be_created");
    symlink(
        &outside,
        f.store_root.join(".locks/refseq.gcf--GCF_000005845.2.lock"),
    )
    .unwrap();
    assert!(f.store().register(&f.source_root, request).is_err());
    assert!(!outside.exists());
    let g = Fixture::new();
    symlink(g.source_root.join("absent"), g.store_root.join("artifacts")).unwrap();
    assert!(g.store().resolve(resolve_refseq("genome_fasta")).is_err());
}

#[test]
fn accessionless_groups_cannot_claim_genome_identity_and_empty_assembly_is_rejected() {
    for empty_canonical in [false, true] {
        let f = Fixture::new();
        let request = f.refseq();
        let root = f.source_root.join("package");
        let path = root.join("ncbi_dataset/data/dataset_catalog.json");
        let mut catalog: serde_json::Value =
            serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        if empty_canonical {
            let moved = catalog["assemblies"][1]["files"].take();
            catalog["assemblies"][1]["files"] = serde_json::json!([]);
            catalog["assemblies"]
                .as_array_mut()
                .unwrap()
                .push(serde_json::json!({"files":moved}));
        } else {
            catalog["assemblies"][1]
                .as_object_mut()
                .unwrap()
                .remove("accession");
            catalog["assemblies"]
                .as_array_mut()
                .unwrap()
                .push(serde_json::json!({"accession":"GCF_000005845.2","files":[]}));
        }
        fs::write(path, serde_json::to_vec_pretty(&catalog).unwrap()).unwrap();
        let error = f.store().register(&f.source_root, request).unwrap_err();
        assert!(matches!(error, StorageError::InvalidPackage(_)));
        assert!(
            error.to_string().contains("accessionless")
                || error
                    .to_string()
                    .contains("canonical assembly catalog group")
        );
        assert!(matches!(f.resolve("genome_fasta"), Resolution::Miss));
    }
}

#[test]
fn raw_json_duplicate_representation_keys_are_rejected() {
    let json = r#"{"source_path":"pdb","requested_ref":"pdb:1CRN","canonical_ref":"pdb:1CRN","declaration":{"provider":"pdb","scope":"entry","representations":{"structure_cif":["first.cif"],"structure_cif":["second.cif"]}}}"#;
    let error = serde_json::from_str::<RegisterRequest>(json).unwrap_err();
    assert!(error.to_string().contains("duplicate representation key"));
    let f = Fixture::new();
    let result = f.store().register(&f.source_root, f.pdb()).unwrap();
    let mut raw = serde_json::to_string(&receipt(&f, &result)).unwrap();
    // Both nested declaration and top-level receipt maps must reject duplicates.
    raw = raw.replace(
        "\"structure_cif\":[\"other.cif\",\"declared.cif\"]",
        "\"structure_cif\":[\"other.cif\"],\"structure_cif\":[\"declared.cif\"]",
    );
    assert!(serde_json::from_str::<Receipt>(&raw)
        .unwrap_err()
        .to_string()
        .contains("duplicate representation key"));
}

#[test]
fn native_backticks_are_legible_in_companion_markdown_paths() {
    let f = Fixture::new();
    let mut request = f.pdb();
    fs::rename(
        f.source_root.join("pdb/declared.cif"),
        f.source_root.join("pdb/weird`name.cif"),
    )
    .unwrap();
    fs::rename(
        f.source_root.join("pdb/other.cif"),
        f.source_root.join("pdb/many``ticks.cif"),
    )
    .unwrap();
    if let NativeDeclaration::Pdb {
        representations, ..
    } = &mut request.declaration
    {
        representations.insert(
            "structure_cif".into(),
            vec!["weird`name.cif".into(), "many``ticks.cif".into()],
        );
    }
    let result = f.store().register(&f.source_root, request).unwrap();
    let readme = fs::read_to_string(
        f.store_root
            .join(result.snapshot.snapshot_path)
            .join("README.md"),
    )
    .unwrap();
    assert!(readme.contains("`` source/weird`name.cif ``"));
    assert!(readme.contains("``` source/many``ticks.cif ```"));
}
#[test]
fn ordinary_genome_example_reads_every_native_fasta_member() {
    let f = Fixture::new();
    let request = f.refseq();
    let root = f.source_root.join("package");
    let catalog_path = root.join("ncbi_dataset/data/dataset_catalog.json");
    fs::write(
        root.join("ncbi_dataset/data/GCF_000005845.2/second.fna"),
        ">second\nAAA\n",
    )
    .unwrap();
    let mut catalog: serde_json::Value =
        serde_json::from_slice(&fs::read(&catalog_path).unwrap()).unwrap();
    catalog["assemblies"][1]["files"].as_array_mut().unwrap().push(serde_json::json!({"filePath":"GCF_000005845.2/second.fna","fileType":"GENOMIC_NUCLEOTIDE_FASTA","uncompressedLengthBytes":"12"}));
    fs::write(catalog_path, serde_json::to_vec_pretty(&catalog).unwrap()).unwrap();
    let inventory = crate::tree::inventory(&root, None).unwrap();
    let mut md5 = String::new();
    for file in inventory
        .iter()
        .filter(|e| e.kind == EntryKind::File && e.path.starts_with("ncbi_dataset/data/"))
    {
        md5.push_str(&format!(
            "{:x}  {}\n",
            Md5::digest(fs::read(root.join(&file.path)).unwrap()),
            file.path
        ));
    }
    fs::write(root.join("md5sum.txt"), md5).unwrap();
    let result = f.store().register(&f.source_root, request).unwrap();
    let readme = fs::read_to_string(
        f.store_root
            .join(&result.snapshot.snapshot_path)
            .join("README.md"),
    )
    .unwrap();
    let example = readme
        .split("```python\n")
        .nth(1)
        .unwrap()
        .split("```")
        .next()
        .unwrap();
    assert!(
        example.contains("genome.fna")
            && example.contains("second.fna")
            && example.contains("for p in map(Path, paths):")
    );
    let Resolution::Ready { paths, .. } = f.resolve("genome_fasta") else {
        panic!("expected ready")
    };
    assert_eq!(paths.len(), 2);
}

#[test]
fn historical_validation_prose_is_allowed_but_changed_method_is_rejected() {
    let f = Fixture::new();
    let request = f.refseq();
    let registered = f.store().register(&f.source_root, request.clone()).unwrap();
    let receipt_path = f.store_root.join(&registered.snapshot.receipt_path);
    let mut historical = receipt(&f, &registered);
    historical.validation.limits =
        vec!["Historical descriptive wording; not an authenticated validation claim".into()];
    fs::write(
        &receipt_path,
        serde_json::to_vec_pretty(&historical).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        f.resolve("genome_fasta"),
        Resolution::Ready { .. }
    ));
    let reused = f.store().register(&f.source_root, request).unwrap();
    assert!(reused.reused);
    assert_eq!(
        receipt(&f, &reused).validation.limits,
        historical.validation.limits
    );
    historical.validation.method = "unperformed_biological_identity_validation".into();
    fs::write(
        receipt_path,
        serde_json::to_vec_pretty(&historical).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        f.resolve("genome_fasta"),
        Resolution::Corrupt { .. }
    ));
}
