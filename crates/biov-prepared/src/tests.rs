use super::*;
use biov_storage::{NativeDeclaration, RegisterRequest};
use md5::{Digest, Md5};
use std::{
    fs,
    sync::{Arc, Barrier},
    thread,
};
use tempfile::TempDir;

const REFERENCE: &str = "refseq.gcf:GCF_000005845.2";
const ACCESSION: &str = "GCF_000005845.2";
const SOURCE_PATH: &str = "ncbi_dataset/data/GCF_000005845.2/genome.fna";
const FASTA: &[u8] = b">0001 description\nACGT\nAC\n>chr2\ttab description\nttNN\n";

struct Fixture {
    _temp: TempDir,
    sources: PathBuf,
    root: PathBuf,
}
impl Fixture {
    fn new() -> Self {
        let temp = tempfile::tempdir().unwrap();
        let sources = temp.path().join("sources");
        let root = temp.path().join("store");
        fs::create_dir(&sources).unwrap();
        fs::create_dir(&root).unwrap();
        Self {
            _temp: temp,
            sources,
            root,
        }
    }
    fn store(&self) -> PreparedStore {
        PreparedStore::new(&self.root).unwrap()
    }
    fn register(&self, fasta: &[u8]) -> PrepareFastaRequest {
        self.register_named(fasta, "genome.fna")
    }
    fn register_named(&self, fasta: &[u8], filename: &str) -> PrepareFastaRequest {
        let package = self.sources.join("package");
        let data = package.join("ncbi_dataset/data");
        let source_path = format!("ncbi_dataset/data/{ACCESSION}/{filename}");
        fs::create_dir_all(data.join(ACCESSION)).unwrap();
        fs::write(
            package.join("README.md"),
            "Native synthetic RefSeq fixture\n",
        )
        .unwrap();
        fs::write(package.join(&source_path), fasta).unwrap();
        fs::write(
            data.join("assembly_data_report.jsonl"),
            format!("{{\"accession\":\"{ACCESSION}\"}}\n"),
        )
        .unwrap();
        let catalog = serde_json::json!({"apiVersion":"V2","assemblies":[
            {"files":[{"filePath":"assembly_data_report.jsonl","fileType":"DATA_REPORT","uncompressedLengthBytes":fs::metadata(data.join("assembly_data_report.jsonl")).unwrap().len().to_string()}]},
            {"accession":ACCESSION,"files":[{"filePath":format!("{ACCESSION}/{filename}"),"fileType":"GENOMIC_NUCLEOTIDE_FASTA","uncompressedLengthBytes":fasta.len().to_string()}]}
        ]});
        fs::write(
            data.join("dataset_catalog.json"),
            serde_json::to_vec_pretty(&catalog).unwrap(),
        )
        .unwrap();
        let mut md5 = String::new();
        for path in [
            "ncbi_dataset/data/dataset_catalog.json".to_owned(),
            "ncbi_dataset/data/assembly_data_report.jsonl".to_owned(),
            source_path.clone(),
        ] {
            md5.push_str(&format!(
                "{:x}  {path}\n",
                Md5::digest(fs::read(package.join(&path)).unwrap())
            ));
        }
        fs::write(package.join("md5sum.txt"), md5).unwrap();
        let registered = NativeStore::new(&self.root)
            .unwrap()
            .register(
                &self.sources,
                RegisterRequest {
                    source_path: "package".into(),
                    requested_ref: REFERENCE.into(),
                    canonical_ref: REFERENCE.into(),
                    declaration: NativeDeclaration::Refseq,
                },
            )
            .unwrap();
        PrepareFastaRequest {
            reference: REFERENCE.into(),
            snapshot_id: registered.snapshot.snapshot_id,
            source_path,
        }
    }
    fn prepare(&self) -> (PrepareFastaRequest, PreparedFasta) {
        let request = self.register(FASTA);
        let result = self.store().prepare_fasta(request.clone()).unwrap();
        (request, result)
    }
}
fn manifest(result: &PreparedFasta) -> FastaProvenance {
    serde_json::from_slice(&fs::read(&result.provenance.execution_host_path).unwrap()).unwrap()
}
fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for member in fs::read_dir(from).unwrap() {
        let member = member.unwrap();
        let target = to.join(member.file_name());
        if member.file_type().unwrap().is_dir() {
            copy_tree(&member.path(), &target);
        } else {
            fs::copy(member.path(), target).unwrap();
        }
    }
}

fn metrics_request(
    fasta: PrepareFastaRequest,
    prepared: &PreparedFasta,
    sequence_id: &str,
    window_size: u64,
) -> FastaWindowMetricsRequest {
    FastaWindowMetricsRequest {
        reference: fasta.reference,
        snapshot_id: fasta.snapshot_id,
        source_path: fasta.source_path,
        recipe_id: prepared.recipe_id.clone(),
        sequence_id: sequence_id.into(),
        window_size,
    }
}

#[test]
fn window_metrics_stream_iupac_case_wrapping_and_same_pass_whole_summary() {
    let fixture = Fixture::new();
    let fasta = b">0001:alt description\r\nacgtry\r\nswkmbd\r\nhvnGCN\r\n>other\r\nGGG";
    let request = fixture.register(fasta);
    let prepared = fixture.store().prepare_fasta(request.clone()).unwrap();
    let selection = fixture
        .store()
        .preflight_fasta_window_metrics(metrics_request(request, &prepared, "0001:alt", 5))
        .unwrap();
    assert_eq!(selection.length, 18);
    assert_eq!(selection.window_row_count, 4);
    assert_eq!(selection.source_bytes, fasta.len() as u64);
    assert_eq!(selection.recipe_id, prepared.recipe_id);
    let mut windows = Vec::new();
    let summary = selection.stream_windows(|row| windows.push(row)).unwrap();
    assert_eq!(
        windows
            .iter()
            .map(|row| (row.start, row.end, row.length, row.is_full_window))
            .collect::<Vec<_>>(),
        vec![
            (0, 5, 5, true),
            (5, 10, 5, true),
            (10, 15, 5, true),
            (15, 18, 3, false)
        ]
    );
    assert_eq!(
        windows
            .iter()
            .map(|row| row.gc_fraction)
            .collect::<Vec<_>>(),
        vec![Some(0.5), None, None, Some(1.0)]
    );
    assert_eq!(
        windows
            .iter()
            .map(|row| row.weighted_gc_fraction)
            .collect::<Vec<_>>(),
        vec![0.5, 0.5, 0.5, 5.0 / 6.0]
    );
    assert!(windows.iter().all(|row| row.sequence_id == "0001:alt"));
    assert_eq!(summary.sequence_id, "0001:alt");
    assert_eq!(summary.length, 18);
    assert_eq!(summary.canonical_base_count, 6);
    assert_eq!(summary.gc_base_count, 4);
    assert_eq!(summary.gc_fraction, Some(2.0 / 3.0));
    assert_eq!(summary.weighted_gc_fraction, 5.0 / 9.0);
    assert_eq!(fs::read(prepared.fasta.execution_host_path).unwrap(), fasta);
}

#[test]
fn window_metrics_exact_multiple_and_all_ambiguous_summary_are_explicit() {
    let fixture = Fixture::new();
    let request = fixture.register(b">allN\nNNnn\nNNnn\n");
    let prepared = fixture.store().prepare_fasta(request.clone()).unwrap();
    let selection = fixture
        .store()
        .preflight_fasta_window_metrics(metrics_request(request, &prepared, "allN", 4))
        .unwrap();
    assert_eq!(selection.window_row_count, 2);
    let mut windows = Vec::new();
    let summary = selection.stream_windows(|row| windows.push(row)).unwrap();
    assert!(windows.iter().all(|row| row.is_full_window));
    assert!(windows.iter().all(|row| row.gc_fraction.is_none()));
    assert_eq!(summary.gc_fraction, None);
    assert_eq!(summary.canonical_base_count, 0);
    assert_eq!(summary.gc_base_count, 0);
    assert_eq!(summary.weighted_gc_fraction, 0.5);
}

#[test]
fn window_metrics_missing_preparation_is_read_only_and_never_created() {
    let fixture = Fixture::new();
    let request = fixture.register(FASTA);
    let store = fixture.store();
    let recipe_id = store.input(&request).unwrap().recipe.id().unwrap();
    let before = fs::read_dir(&fixture.root).unwrap().count();
    let result = store.preflight_fasta_window_metrics(FastaWindowMetricsRequest {
        reference: request.reference,
        snapshot_id: request.snapshot_id,
        source_path: request.source_path,
        recipe_id,
        sequence_id: "0001".into(),
        window_size: 2,
    });
    assert!(matches!(result, Err(PreparedError::InputUnavailable(_))));
    assert_eq!(fs::read_dir(&fixture.root).unwrap().count(), before);
    assert!(!fixture.root.join("prepared").exists());
    // Native registration already created its own empty staging parent.
    assert_eq!(
        fs::read_dir(fixture.root.join(".staging")).unwrap().count(),
        0
    );
}

#[test]
fn window_metrics_request_deserialization_is_flat_and_rejects_unknown_fields() {
    let value = serde_json::json!({
        "reference": REFERENCE, "snapshot_id": format!("sha256-{}", "a".repeat(64)),
        "source_path": SOURCE_PATH, "recipe_id": format!("sha256-{}", "b".repeat(64)),
        "sequence_id": "0001", "window_size": 3,
    });
    let decoded: FastaWindowMetricsRequest = serde_json::from_value(value.clone()).unwrap();
    assert_eq!(decoded.reference, REFERENCE);
    assert_eq!(decoded.window_size, 3);
    assert_eq!(serde_json::to_value(decoded).unwrap(), value);
    let mut unknown = value;
    unknown["unexpected"] = serde_json::json!(true);
    assert!(serde_json::from_value::<FastaWindowMetricsRequest>(unknown).is_err());
}

#[test]
fn window_metrics_invalid_identity_selection_and_bounds_reject() {
    let fixture = Fixture::new();
    let (request, prepared) = fixture.prepare();
    let valid = metrics_request(request, &prepared, "0001", 2);
    for recipe in [
        "sha256-zero".to_owned(),
        format!("sha256-{}", "0".repeat(64)),
        format!("sha256-{}", "A".repeat(64)),
    ] {
        let mut request = valid.clone();
        request.recipe_id = recipe;
        assert!(matches!(
            fixture.store().preflight_fasta_window_metrics(request),
            Err(PreparedError::InvalidInput(_))
        ));
    }
    for id in [
        "",
        "0001 description",
        "0001:1-2",
        "00001",
        "Chr2",
        "chr2\n",
        "é",
    ] {
        let mut request = valid.clone();
        request.sequence_id = id.into();
        assert!(matches!(
            fixture.store().preflight_fasta_window_metrics(request),
            Err(PreparedError::InvalidInput(_))
        ));
    }
    let mut request = valid.clone();
    request.reference = "refseq.gcf:GCF_000005845".into();
    assert!(matches!(
        fixture.store().preflight_fasta_window_metrics(request),
        Err(PreparedError::InvalidInput(_))
    ));
    let mut request = valid.clone();
    request.snapshot_id = "latest".into();
    assert!(matches!(
        fixture.store().preflight_fasta_window_metrics(request),
        Err(PreparedError::InvalidInput(_))
    ));
    let mut request = valid.clone();
    request.source_path = "../outside.fna".into();
    assert!(matches!(
        fixture.store().preflight_fasta_window_metrics(request),
        Err(PreparedError::InvalidInput(_))
    ));
    let mut request = valid.clone();
    request.source_path = "ncbi_dataset/data/assembly_data_report.jsonl".into();
    assert!(matches!(
        fixture.store().preflight_fasta_window_metrics(request),
        Err(PreparedError::InvalidInput(_))
    ));
    let mut request = valid.clone();
    request.window_size = 0;
    assert!(matches!(
        fixture.store().preflight_fasta_window_metrics(request),
        Err(PreparedError::InvalidInput(_))
    ));
    let mut request = valid;
    request.window_size = MAX_WINDOW_BASES + 1;
    assert!(matches!(
        fixture.store().preflight_fasta_window_metrics(request),
        Err(PreparedError::Limit(_))
    ));
}

#[test]
fn window_metrics_source_and_prepared_corruption_reject_before_and_after_emission() {
    for corrupt_source in [false, true] {
        let fixture = Fixture::new();
        let (request, prepared) = fixture.prepare();
        let selection = fixture
            .store()
            .preflight_fasta_window_metrics(metrics_request(request.clone(), &prepared, "0001", 2))
            .unwrap();
        let mut emitted = 0;
        let result = selection.stream_windows(|_| {
            emitted += 1;
            if emitted == 1 {
                if corrupt_source {
                    let mut source = FASTA.to_vec();
                    source[18] = b'T';
                    fs::write(&prepared.fasta.execution_host_path, source).unwrap();
                } else {
                    fs::write(&prepared.dictionary.execution_host_path, b"modified").unwrap();
                }
            }
        });
        assert!(
            matches!(result, Err(PreparedError::Corrupt(_))),
            "{result:?}"
        );
        assert!(fixture
            .store()
            .preflight_fasta_window_metrics(metrics_request(request, &prepared, "0001", 2),)
            .is_err());
    }
}

#[test]
fn window_metrics_caps_each_query_across_a_sequence_larger_than_the_cap() {
    let fixture = Fixture::new();
    let mut fasta = b">long\n".to_vec();
    let line = vec![b'a'; IO_BUFFER_BYTES];
    for _ in 0..MAX_WINDOW_BASES as usize / IO_BUFFER_BYTES + 1 {
        fasta.extend_from_slice(&line);
        fasta.push(b'\n');
    }
    fasta.extend_from_slice(b"gcn");
    let request = fixture.register(&fasta);
    let prepared = fixture.store().prepare_fasta(request.clone()).unwrap();
    let selection = fixture
        .store()
        .preflight_fasta_window_metrics(metrics_request(
            request,
            &prepared,
            "long",
            MAX_WINDOW_BASES,
        ))
        .unwrap();
    assert_eq!(selection.window_row_count, 2);
    let mut windows = Vec::new();
    let summary = selection.stream_windows(|row| windows.push(row)).unwrap();
    assert_eq!(windows[0].length, MAX_WINDOW_BASES);
    assert!(windows[0].is_full_window);
    assert_eq!(windows[1].length, IO_BUFFER_BYTES as u64 + 3);
    assert!(!windows[1].is_full_window);
    assert_eq!(
        summary.length,
        MAX_WINDOW_BASES + IO_BUFFER_BYTES as u64 + 3
    );
    assert_eq!(summary.gc_base_count, 2);
    assert_eq!(summary.canonical_base_count, summary.length - 1);
    assert_eq!(
        summary.weighted_gc_fraction,
        15.0 / (6.0 * summary.length as f64)
    );
}
#[test]
fn construction_is_read_only_and_requires_existing_trusted_root() {
    let fixture = Fixture::new();
    fixture.store();
    assert_eq!(fs::read_dir(&fixture.root).unwrap().count(), 0);
    assert!(PreparedStore::new(fixture.root.join("missing")).is_err());
}
#[test]
fn prepare_preserves_native_bytes_publishes_conventional_outputs_and_reuses() {
    let fixture = Fixture::new();
    let (request, result) = fixture.prepare();
    assert!(!result.reused);
    assert_eq!(result.sequence_count, 2);
    assert_eq!(result.total_bases, 10);
    assert_eq!(fs::read(&result.fasta.execution_host_path).unwrap(), FASTA);
    assert_eq!(
        fs::read(fixture.sources.join("package").join(SOURCE_PATH)).unwrap(),
        FASTA
    );
    assert_eq!(
        fs::read(&result.dictionary.execution_host_path).unwrap(),
        b"sequence_id\tlength\n0001\t6\nchr2\t4\n"
    );
    let parsed = noodles_fasta::fai::io::Reader::new(
        fsutil::read_bounded(
            Path::new(&result.fai.execution_host_path),
            MAX_METADATA_BYTES,
        )
        .unwrap()
        .as_slice(),
    )
    .read_index()
    .unwrap();
    assert_eq!(parsed.as_ref()[0].name(), b"0001");
    assert_eq!(parsed.as_ref()[0].length(), 6);
    let saved = manifest(&result);
    assert_eq!(saved.recipe_id, saved.recipe.id().unwrap());
    assert_eq!(
        saved.recipe.source_sha256,
        format!("{:x}", Sha256::digest(FASTA))
    );
    assert_eq!(saved.outputs.len(), 3);
    let before = fs::read(&result.provenance.execution_host_path).unwrap();
    let reused = fixture.store().prepare_fasta(request).unwrap();
    assert!(reused.reused);
    assert_eq!(reused.recipe_id, result.recipe_id);
    assert_eq!(
        before,
        fs::read(reused.provenance.execution_host_path).unwrap()
    );
    assert_eq!(
        fs::read_dir(fixture.root.join(".staging")).unwrap().count(),
        0
    );
}
#[test]
fn request_requires_exact_canonical_version_pin_and_selected_representation_path() {
    let fixture = Fixture::new();
    let request = fixture.register(FASTA);
    for reference in [
        "refseq.gcf:GCF_000005845",
        "GCF_000005845.2",
        "refseq.gcf://GCF_000005845.2",
        "REFSEQ.GCF:GCF_000005845.2",
        "uniprot:P12345",
    ] {
        let mut invalid = request.clone();
        invalid.reference = reference.into();
        assert!(matches!(
            fixture.store().prepare_fasta(invalid),
            Err(PreparedError::InvalidInput(_))
        ));
    }
    for pin in [
        "latest",
        "sha256-abcd",
        "sha256-AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    ] {
        let mut invalid = request.clone();
        invalid.snapshot_id = pin.into();
        assert!(matches!(
            fixture.store().prepare_fasta(invalid),
            Err(PreparedError::InvalidInput(_))
        ));
    }
    for path in [
        "../genome.fna",
        "/tmp/genome.fna",
        "./genome.fna",
        "a//b",
        "a\\b",
        "ncbi_dataset/data/dataset_catalog.json",
        "README.md",
        "NUL",
        "a/../b",
    ] {
        let mut invalid = request.clone();
        invalid.source_path = path.into();
        assert!(
            matches!(
                fixture.store().prepare_fasta(invalid),
                Err(PreparedError::InvalidInput(_))
            ),
            "{path}"
        );
    }
    assert!(!fixture.root.join("prepared").exists());
}
#[test]
fn missing_snapshot_is_explicit_without_creating_prepared_directories() {
    let fixture = Fixture::new();
    let request = PrepareFastaRequest {
        reference: REFERENCE.into(),
        snapshot_id: format!("sha256-{}", "0".repeat(64)),
        source_path: SOURCE_PATH.into(),
    };
    assert!(matches!(
        fixture.store().prepare_fasta(request),
        Err(PreparedError::InputUnavailable(_))
    ));
    assert_eq!(fs::read_dir(&fixture.root).unwrap().count(), 0);
}
#[test]
fn malformed_compressed_and_unsupported_fasta_never_publish() {
    for source in [
        b"\x1f\x8bgarbage".as_slice(),
        b">a\nAC\n>a\nGT\n",
        b">a\nA C\n",
        b">a\nA\nAAAA\n",
        b">a\n",
        b"not fasta\n",
    ] {
        let fixture = Fixture::new();
        let request = fixture.register(source);
        assert!(matches!(
            fixture.store().prepare_fasta(request),
            Err(PreparedError::InvalidFasta(_))
        ));
        assert_eq!(
            fs::read_dir(fixture.root.join("prepared")).unwrap().count(),
            0
        );
        assert_eq!(
            fs::read_dir(fixture.root.join(".staging")).unwrap().count(),
            0
        );
    }
    let fixture = Fixture::new();
    let request = fixture.register_named(FASTA, "genome.fna.gz");
    assert!(matches!(
        fixture.store().prepare_fasta(request),
        Err(PreparedError::InvalidFasta(_))
    ));
}
#[test]
fn changed_or_missing_native_input_never_reuses_existing_output() {
    for remove in [false, true] {
        let fixture = Fixture::new();
        let (request, result) = fixture.prepare();
        let before = fs::read(&result.provenance.execution_host_path).unwrap();
        if remove {
            fs::remove_file(&result.fasta.execution_host_path).unwrap();
        } else {
            fs::write(&result.fasta.execution_host_path, b">mutated\nAAAA\n").unwrap();
        }
        assert!(fixture.store().prepare_fasta(request).is_err());
        assert_eq!(
            before,
            fs::read(&result.provenance.execution_host_path).unwrap()
        );
    }
}
#[test]
fn every_output_tamper_missing_member_and_extra_member_is_rejected_without_replacement() {
    for choice in 0..9 {
        let fixture = Fixture::new();
        let (request, result) = fixture.prepare();
        let dir = Path::new(&result.provenance.execution_host_path)
            .parent()
            .unwrap();
        let path = match choice {
            0 | 4 => &result.fai.execution_host_path,
            1 | 5 => &result.dictionary.execution_host_path,
            2 | 6 => &result.readme.execution_host_path,
            _ => &result.provenance.execution_host_path,
        };
        if choice < 4 {
            fs::write(path, "tampered\n").unwrap();
        } else if choice < 8 {
            fs::remove_file(path).unwrap();
        } else {
            fs::write(dir.join("unexpected.txt"), "extra").unwrap();
        }
        assert!(
            fixture.store().prepare_fasta(request).is_err(),
            "tamper choice {choice}"
        );
    }
}
#[test]
fn tampering_both_index_and_its_provenance_hash_cannot_forge_reuse() {
    let fixture = Fixture::new();
    let (request, result) = fixture.prepare();
    let changed = b"0001\t999\t18\t4\t5\nchr2\t4\t49\t4\t5\n";
    fs::write(&result.fai.execution_host_path, changed).unwrap();
    let mut record = manifest(&result);
    record.outputs[0].bytes = changed.len() as u64;
    record.outputs[0].sha256 = format!("{:x}", Sha256::digest(changed));
    fs::write(
        &result.provenance.execution_host_path,
        serde_json::to_vec_pretty(&record).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        fixture.store().prepare_fasta(request),
        Err(PreparedError::Corrupt(_))
    ));
    assert_eq!(fs::read(&result.fai.execution_host_path).unwrap(), changed);
}
#[test]
fn recipe_identity_changes_for_inputs_tools_versions_and_parameters_not_output_hashes() {
    let fixture = Fixture::new();
    let (_, result) = fixture.prepare();
    let original = manifest(&result).recipe;
    let mut variants = Vec::new();
    let mut altered = original.clone();
    altered.algorithm_revision += 1;
    variants.push(altered);
    let mut altered = original.clone();
    altered.implementation_version = "0.2.0".into();
    variants.push(altered);
    let mut altered = original.clone();
    altered.library_version = "0.67.0".into();
    variants.push(altered);
    let mut altered = original.clone();
    altered.source_sha256 = "0".repeat(64);
    variants.push(altered);
    let mut altered = original.clone();
    altered.source_path = "other/genome.fna".into();
    variants.push(altered);
    let mut altered = original.clone();
    altered.snapshot_id = format!("sha256-{}", "0".repeat(64));
    variants.push(altered);
    let mut altered = original.clone();
    altered.parameters.max_physical_line_bytes += 1;
    variants.push(altered);
    let mut altered = original.clone();
    altered.parameters.dictionary_schema.push_str("-changed");
    variants.push(altered);
    for altered in variants {
        assert_ne!(original.id().unwrap(), altered.id().unwrap());
    }
    let mut record = manifest(&result);
    record.outputs[0].sha256 = "0".repeat(64);
    assert_eq!(original.id().unwrap(), record.recipe.id().unwrap());
}
#[test]
fn changing_native_package_bytes_creates_new_recipe_without_replacing_old() {
    let fixture = Fixture::new();
    let (_, first) = fixture.prepare();
    // A new source-byte snapshot, even when the biological accession is unchanged.
    let request = fixture.register(b">0001\nACGT\nACGT\n");
    let second = fixture.store().prepare_fasta(request).unwrap();
    assert_ne!(first.snapshot_id, second.snapshot_id);
    assert_ne!(first.recipe_id, second.recipe_id);
    assert!(Path::new(&first.fai.execution_host_path).is_file());
    assert_eq!(fs::read(&first.fasta.execution_host_path).unwrap(), FASTA);
}
#[test]
fn exact_portable_closure_can_move_without_original_source_or_store() {
    let fixture = Fixture::new();
    let (request, result) = fixture.prepare();
    let moved = fixture._temp.path().join("unrelated-reader-root");
    fs::create_dir(&moved).unwrap();
    let p = manifest(&result);
    let snapshot = p
        .source_snapshot_relative_path
        .strip_prefix("../../")
        .unwrap();
    copy_tree(&fixture.root.join(snapshot), &moved.join(snapshot));
    let prepared = format!("prepared/{}", result.recipe_id);
    copy_tree(&fixture.root.join(&prepared), &moved.join(prepared));
    fs::remove_dir_all(&fixture.sources).unwrap();
    fs::remove_dir_all(&fixture.root).unwrap();
    let moved_store = PreparedStore::new(&moved).unwrap();
    let before = fs::read_dir(&moved).unwrap().count();
    let selection = moved_store
        .preflight_fasta_window_metrics(metrics_request(request.clone(), &result, "0001", 4))
        .unwrap();
    let mut windows = Vec::new();
    let summary = selection.stream_windows(|row| windows.push(row)).unwrap();
    assert_eq!(summary.length, 6);
    assert_eq!(summary.gc_fraction, Some(0.5));
    assert_eq!(summary.weighted_gc_fraction, 0.5);
    assert_eq!(windows.len(), 2);
    assert_eq!((windows[1].start, windows[1].end), (4, 6));
    assert_eq!(fs::read_dir(&moved).unwrap().count(), before);
    let moved_result = moved_store.prepare_fasta(request).unwrap();
    assert!(moved_result.reused);
    assert_eq!(moved_result.recipe_id, result.recipe_id);
    assert_eq!(
        fs::read(&moved_result.fasta.execution_host_path).unwrap(),
        FASTA
    );
    let text = fs::read_to_string(&moved_result.provenance.execution_host_path).unwrap();
    assert!(!text.contains(fixture._temp.path().to_str().unwrap()));
}
#[test]
fn source_paths_with_quotes_backticks_unicode_and_spaces_remain_safe_and_exact() {
    let fixture = Fixture::new();
    let request = fixture.register_named(FASTA, "genome 'quote` ü space.fna");
    let result = fixture.store().prepare_fasta(request).unwrap();
    let p = manifest(&result);
    assert!(p
        .source_relative_path
        .ends_with("genome 'quote` ü space.fna"));
    let readme = fs::read_to_string(&result.readme.execution_host_path).unwrap();
    assert!(!readme.contains("genome 'quote` ü space.fna"));
    assert!(readme.contains("source = Path(p[\"source_relative_path\"])"));
}
#[test]
fn stale_interrupted_stage_is_invisible_and_preserved_until_explicit_cleanup() {
    let fixture = Fixture::new();
    let request = fixture.register(FASTA);
    let interrupted = fixture.root.join(".staging/prepare-fasta-interrupted");
    fs::create_dir_all(&interrupted).unwrap();
    fs::write(interrupted.join("sequences.fai"), "partial").unwrap();
    let result = fixture.store().prepare_fasta(request).unwrap();
    assert!(!result.reused);
    assert_eq!(
        fs::read(interrupted.join("sequences.fai")).unwrap(),
        b"partial"
    );
    assert_eq!(
        fs::read_dir(fixture.root.join("prepared")).unwrap().count(),
        1
    );
}
#[test]
fn concurrent_identical_requests_publish_one_verified_winner() {
    let fixture = Fixture::new();
    let request = fixture.register(FASTA);
    let barrier = Arc::new(Barrier::new(4));
    let mut workers = Vec::new();
    for _ in 0..4 {
        let barrier = barrier.clone();
        let request = request.clone();
        let store = fixture.store();
        workers.push(thread::spawn(move || {
            barrier.wait();
            store.prepare_fasta(request).unwrap()
        }));
    }
    let results: Vec<_> = workers
        .into_iter()
        .map(|worker| worker.join().unwrap())
        .collect();
    assert_eq!(results.iter().filter(|result| !result.reused).count(), 1);
    assert!(results
        .iter()
        .all(|result| result.recipe_id == results[0].recipe_id));
    assert_eq!(
        fs::read_dir(fixture.root.join("prepared")).unwrap().count(),
        1
    );
}
#[cfg(unix)]
#[test]
fn links_at_root_prepared_target_outputs_and_locks_are_rejected() {
    use std::os::unix::fs::symlink;
    let fixture = Fixture::new();
    let alias = fixture._temp.path().join("alias");
    symlink(&fixture.root, &alias).unwrap();
    assert!(PreparedStore::new(alias).is_err());
    for location in ["prepared", ".staging", ".locks"] {
        let fixture = Fixture::new();
        let request = fixture.register(FASTA);
        let path = fixture.root.join(location);
        if path.exists() {
            fs::remove_dir_all(&path).unwrap();
        }
        symlink(&fixture.sources, path).unwrap();
        assert!(fixture.store().prepare_fasta(request).is_err());
    }
    for name in [
        "sequences.fai",
        "sequences.tsv",
        "provenance.json",
        "README.md",
    ] {
        let fixture = Fixture::new();
        let (request, result) = fixture.prepare();
        let path = Path::new(&result.provenance.execution_host_path)
            .parent()
            .unwrap()
            .join(name);
        let outside = fixture._temp.path().join("outside");
        fs::copy(&path, &outside).unwrap();
        fs::remove_file(&path).unwrap();
        symlink(outside, path).unwrap();
        assert!(fixture.store().prepare_fasta(request).is_err());
    }
}
#[test]
fn no_replace_publication_preserves_an_existing_target() {
    let root = tempfile::tempdir().unwrap();
    let source = root.path().join("source");
    let target = root.path().join("target");
    fs::create_dir(&source).unwrap();
    fs::create_dir(&target).unwrap();
    fs::write(target.join("keep"), "winner").unwrap();
    assert!(fsutil::publish(&source, &target).is_err());
    assert_eq!(fs::read(target.join("keep")).unwrap(), b"winner");
    assert!(source.is_dir());
}
#[test]
fn oversized_provenance_is_rejected_before_deserialization() {
    let fixture = Fixture::new();
    let (request, result) = fixture.prepare();
    fs::write(
        &result.provenance.execution_host_path,
        vec![b' '; MAX_PROVENANCE_BYTES + 1],
    )
    .unwrap();
    assert!(matches!(
        fixture.store().prepare_fasta(request),
        Err(PreparedError::Limit(_))
    ));
}

#[cfg(unix)]
#[test]
fn oversized_host_path_response_fails_before_publication() {
    let fixture = Fixture::new();
    let request = fixture.register(FASTA);
    let mut long_parent = fixture._temp.path().join("long-root");
    for _ in 0..10 {
        long_parent.push("\u{1}".repeat(220));
    }
    fs::create_dir_all(&long_parent).unwrap();
    let root = long_parent.join("store");
    fs::rename(&fixture.root, &root).unwrap();
    let store = PreparedStore::new(&root).unwrap();
    assert!(matches!(
        store.prepare_fasta(request),
        Err(PreparedError::Limit(_))
    ));
    assert_eq!(fs::read_dir(root.join("prepared")).unwrap().count(), 0);
    assert_eq!(fs::read_dir(root.join(".staging")).unwrap().count(), 0);
}

#[test]
fn ordinary_source_edit_during_indexing_prevents_publication() {
    use std::{
        io::Write,
        time::{Duration, Instant},
    };
    let fixture = Fixture::new();
    let mut fasta = b">chr\n".to_vec();
    let line = [b"A".repeat(80), b"\n".to_vec()].concat();
    for _ in 0..400_000 {
        fasta.extend_from_slice(&line);
    }
    let request = fixture.register(&fasta);
    let snapshot_path = fixture.root.join(format!(
        "artifacts/refseq.gcf/{ACCESSION}/snapshots/{}/source/{SOURCE_PATH}",
        request.snapshot_id
    ));
    let store = fixture.store();
    let worker = thread::spawn(move || store.prepare_fasta(request));
    let deadline = Instant::now() + Duration::from_secs(30);
    loop {
        let staged = fs::read_dir(fixture.root.join(".staging"))
            .unwrap()
            .any(|entry| {
                entry
                    .unwrap()
                    .file_name()
                    .to_string_lossy()
                    .starts_with("prepare-fasta-")
            });
        if staged {
            break;
        }
        assert!(
            !worker.is_finished(),
            "preparation finished before observing its indexing stage"
        );
        assert!(Instant::now() < deadline, "indexing stage did not appear");
        thread::sleep(Duration::from_millis(1));
    }
    // Append valid sequence bytes. Whether the open stream has passed EOF or not,
    // its digest or final complete native revalidation must detect this edit.
    fs::OpenOptions::new()
        .append(true)
        .open(snapshot_path)
        .unwrap()
        .write_all(b"A\n")
        .unwrap();
    assert!(worker.join().unwrap().is_err());
    assert_eq!(
        fs::read_dir(fixture.root.join("prepared")).unwrap().count(),
        0
    );
}
