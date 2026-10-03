//! Durable native-export reuse, independent of old session handles.
use biov_data::{DatasetStore, ExportRequest, OpenRequest, ReopenRequest};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{fs, path::Path};

fn fixture(csv: &str, schema: Value) -> (tempfile::TempDir, Value, Value) {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("table.csv"), csv).unwrap();
    let mut store = DatasetStore::new(dir.path(), dir.path()).unwrap();
    let opened = store.open(serde_json::from_value::<OpenRequest>(json!({
        "path":"table.csv", "schema":schema,
        "metadata":{"identifier":"uniprot:P05067", "reference":"GRCh38", "coordinates":"1-based closed", "units":"reads"}
    })).unwrap()).unwrap();
    let exported = store
        .export(ExportRequest {
            dataset_id: opened["dataset_id"].as_str().unwrap().into(),
        })
        .unwrap();
    (dir, opened, exported)
}
fn name(export: &Value) -> String {
    Path::new(export["record_path"].as_str().unwrap())
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .into()
}
fn reopen(dir: &Path, record: &str) -> biov_data::Result<Value> {
    DatasetStore::new(dir, dir).unwrap().reopen(ReopenRequest {
        record_path: record.into(),
        preview_rows: 5,
    })
}
fn rewrite_record(export: &Value, record: &Value) {
    fs::write(
        export["record_path"].as_str().unwrap(),
        serde_json::to_vec(record).unwrap(),
    )
    .unwrap();
}
#[test]
fn typed_nulls_empty_strings_long_text_and_large_integers_survive_restart() {
    let long = "基因".repeat(200);
    let csv = format!("id,count,ratio,passed,description\n00123,9007199254740993,1.25,true,{long}\n00456,-9223372036854775808,-2.5,false,\"\"\n00000,,,,\n");
    let (dir, opened, exported) = fixture(
        &csv,
        json!({"count":"int64","ratio":"float64","passed":"boolean"}),
    );
    let mut restarted = DatasetStore::new(dir.path(), dir.path()).unwrap();
    assert!(restarted
        .preview(serde_json::from_value(json!({"dataset_id":opened["dataset_id"]})).unwrap())
        .is_err());
    let recovered = restarted
        .reopen(ReopenRequest {
            record_path: name(&exported),
            preview_rows: 5,
        })
        .unwrap();
    assert_ne!(opened["dataset_id"], recovered["dataset_id"]);
    for field in [
        "row_count",
        "schema",
        "preview",
        "scientific_metadata",
        "source_sha256",
    ] {
        assert_eq!(opened[field], recovered[field], "field {field}");
    }
    assert_eq!(recovered["scientific_metadata"]["species"], Value::Null);
    assert_eq!(recovered["preview"]["rows"][1][4], "");
    assert_eq!(
        recovered["preview"]["rows"][2],
        json!(["00000", null, null, null, null])
    );
    let verification = &recovered["reopen_verification"];
    assert_eq!(verification["artifact_sha256"], exported["sha256"]);
    assert_eq!(
        verification["original_provenance"],
        "recorded_claims_not_independently_verified"
    );
    assert_eq!(verification["authenticity"], "not_established");
    assert_eq!(
        verification["record_sha256"],
        format!(
            "{:x}",
            Sha256::digest(fs::read(exported["record_path"].as_str().unwrap()).unwrap())
        )
    );
    // Read a nonpreview cell through a complete query to confirm no truncation.
    let query = restarted.query(serde_json::from_value(json!({"dataset_id":recovered["dataset_id"],"filter":{"column":"description","op":"eq","value":long},"select":["id"],"preview_rows":5})).unwrap()).unwrap();
    assert_eq!(query["row_count"], 1);
    assert_eq!(query["preview"]["rows"], json!([["00123"]]));
    let reexported = restarted
        .export(ExportRequest {
            dataset_id: query["dataset_id"].as_str().unwrap().into(),
        })
        .unwrap();
    assert_eq!(reexported["record"]["reopen_verification"], *verification);
    let rerecovered = reopen(dir.path(), &name(&reexported)).unwrap();
    assert_eq!(rerecovered["preview"]["rows"], query["preview"]["rows"]);
}
#[test]
fn original_version_one_record_and_empty_typed_exports_reopen() {
    let (dir, _, exported) = fixture("id,count\n001,1\n", json!({"count":"int64"}));
    let mut old = exported["record"].clone();
    old["record_version"] = json!(1);
    old.as_object_mut().unwrap().remove("reopen_verification");
    rewrite_record(&exported, &old);
    let opened = reopen(dir.path(), &name(&exported)).unwrap();
    assert_eq!(opened["reopen_verification"]["record_version"], 1);
    let mut store = DatasetStore::new(dir.path(), dir.path()).unwrap();
    let source = store
        .reopen(ReopenRequest {
            record_path: name(&exported),
            preview_rows: 0,
        })
        .unwrap();
    let empty = store.query(serde_json::from_value(json!({"dataset_id":source["dataset_id"],"filter":{"column":"count","op":"gt","value":10}})).unwrap()).unwrap();
    let empty_export = store
        .export(ExportRequest {
            dataset_id: empty["dataset_id"].as_str().unwrap().into(),
        })
        .unwrap();
    let recovered = reopen(dir.path(), &name(&empty_export)).unwrap();
    assert_eq!(recovered["row_count"], 0);
    assert_eq!(recovered["schema"], source["schema"]);
}
#[test]
fn repeated_reopen_reexport_preserves_bounded_provenance_without_nesting() {
    let (dir, _, mut exported) = fixture("id,value\n001,1\n002,2\n", json!({}));
    let original = exported["record"]["provenance"].clone();
    for _ in 0..20 {
        let mut store = DatasetStore::new(dir.path(), dir.path()).unwrap();
        let source = store
            .reopen(ReopenRequest {
                record_path: name(&exported),
                preview_rows: 1,
            })
            .unwrap();
        exported = store
            .export(ExportRequest {
                dataset_id: source["dataset_id"].as_str().unwrap().into(),
            })
            .unwrap();
        assert_eq!(exported["record"]["provenance"], original);
        assert!(exported["record"]["reopen_verification"]
            .get("reopen_verification")
            .is_none());
        assert!(serde_json::to_vec(&exported["record"]).unwrap().len() < 4096);
    }
}
#[test]
fn strict_record_schema_metadata_and_claim_validation() {
    let (dir, _, export) = fixture("id,count\n001,1\n", json!({"count":"int64"}));
    let original = export["record"].clone();
    for (pointer, value) in [
        ("/record_version", json!(77)),
        ("/format", json!("parquet")),
        ("/bytes", json!(1)),
        ("/sha256", json!("f".repeat(64))),
        ("/row_count", json!(2)),
        ("/schema/0/dtype", json!("i64")),
        ("/schema/0/name", json!("wrong")),
        (
            "/scientific_metadata/identifier",
            json!("not an identifier"),
        ),
        ("/scientific_metadata/species", json!("x".repeat(257))),
        ("/metadata_status", json!("verified")),
        ("/provenance/source_sha256", json!("invalid")),
        ("/provenance/input_consistency", json!("verified")),
        ("/software/polars", json!("unknown")),
    ] {
        let mut record = original.clone();
        *record.pointer_mut(pointer).unwrap() = value;
        rewrite_record(&export, &record);
        assert!(
            reopen(dir.path(), &name(&export)).is_err(),
            "accepted bad {pointer}"
        );
    }
    for pointer in [
        "",
        "/scientific_metadata",
        "/provenance",
        "/software",
        "/schema/0",
    ] {
        let mut record = original.clone();
        record
            .pointer_mut(pointer)
            .unwrap()
            .as_object_mut()
            .unwrap()
            .insert("unexpected".into(), json!(true));
        rewrite_record(&export, &record);
        assert!(
            reopen(dir.path(), &name(&export)).is_err(),
            "accepted unknown key at {pointer}"
        );
    }
    for field in [
        "record_version",
        "format",
        "bytes",
        "sha256",
        "row_count",
        "schema",
        "scientific_metadata",
        "metadata_status",
        "provenance",
        "software",
    ] {
        let mut record = original.clone();
        record.as_object_mut().unwrap().remove(field);
        rewrite_record(&export, &record);
        assert!(
            reopen(dir.path(), &name(&export)).is_err(),
            "accepted missing {field}"
        );
    }
    fs::write(
        export["record_path"].as_str().unwrap(),
        "{\"record_version\":2,\"record_version\":1}",
    )
    .unwrap();
    assert!(reopen(dir.path(), &name(&export)).is_err());
}
#[test]
fn missing_files_tampering_and_resource_caps_fail_explicitly() {
    let (dir, _, export) = fixture("id\n001\n", json!({}));
    let record_name = name(&export);
    assert!(reopen(dir.path(), "missing.json").is_err());
    assert!(reopen(dir.path(), "table.csv").is_err());
    let original = export["record"].clone();
    let mut record = original.clone();
    record["row_count"] = json!(usize::MAX);
    rewrite_record(&export, &record);
    assert!(reopen(dir.path(), &record_name)
        .unwrap_err()
        .to_string()
        .contains("allocation"));
    fs::File::create(export["record_path"].as_str().unwrap())
        .unwrap()
        .set_len(64 * 1024 + 1)
        .unwrap();
    assert!(reopen(dir.path(), &record_name)
        .unwrap_err()
        .to_string()
        .contains("byte limit"));
    rewrite_record(&export, &original);
    let ipc = export["execution_host_path"].as_str().unwrap();
    let original_bytes = fs::read(ipc).unwrap();
    fs::write(ipc, b"changed").unwrap();
    assert!(reopen(dir.path(), &record_name)
        .unwrap_err()
        .to_string()
        .contains("size"));
    let mut changed = original_bytes.clone();
    changed[10] ^= 1;
    fs::write(ipc, &changed).unwrap();
    assert!(reopen(dir.path(), &record_name)
        .unwrap_err()
        .to_string()
        .contains("SHA-256"));
    fs::File::create(ipc)
        .unwrap()
        .set_len(65 * 1024 * 1024 + 1)
        .unwrap();
    assert!(reopen(dir.path(), &record_name)
        .unwrap_err()
        .to_string()
        .contains("byte limit"));
    fs::remove_file(ipc).unwrap();
    assert!(reopen(dir.path(), &record_name).is_err());
    fs::write(ipc, original_bytes).unwrap();
    fs::remove_file(export["record_path"].as_str().unwrap()).unwrap();
    assert!(reopen(dir.path(), &record_name).is_err());
}
#[test]
fn data_root_confinement_and_same_directory_ipc_are_enforced() {
    let (dir, _, export) = fixture("id\n001\n", json!({}));
    for path in ["../outside.json", "/absolute.json", "a/../record.json", ""] {
        assert!(reopen(dir.path(), path).is_err());
    }
    fs::create_dir(dir.path().join("nested")).unwrap();
    let record_name = name(&export);
    fs::rename(
        export["record_path"].as_str().unwrap(),
        dir.path().join("nested").join(&record_name),
    )
    .unwrap();
    fs::rename(
        export["execution_host_path"].as_str().unwrap(),
        dir.path()
            .join("nested")
            .join(export["record"]["file"].as_str().unwrap()),
    )
    .unwrap();
    assert!(reopen(dir.path(), &format!("nested/{record_name}")).is_ok());
    let record_path = dir.path().join("nested").join(&record_name);
    for bad in ["../a.arrow", "/a.arrow", "subdir/a.arrow"] {
        let mut record = export["record"].clone();
        record["file"] = json!(bad);
        fs::write(&record_path, serde_json::to_vec(&record).unwrap()).unwrap();
        assert!(reopen(dir.path(), &format!("nested/{record_name}")).is_err());
    }
}
#[cfg(unix)]
#[test]
fn symlink_escapes_for_record_and_ipc_are_rejected() {
    use std::os::unix::fs::symlink;
    let (outside, _, export) = fixture("id\n001\n", json!({}));
    let inside = tempfile::tempdir().unwrap();
    symlink(
        export["record_path"].as_str().unwrap(),
        inside.path().join(name(&export)),
    )
    .unwrap();
    assert!(reopen(inside.path(), &name(&export))
        .unwrap_err()
        .to_string()
        .contains("outside"));
    fs::remove_file(inside.path().join(name(&export))).unwrap();
    fs::copy(
        export["record_path"].as_str().unwrap(),
        inside.path().join(name(&export)),
    )
    .unwrap();
    symlink(
        outside
            .path()
            .join(export["record"]["file"].as_str().unwrap()),
        inside
            .path()
            .join(export["record"]["file"].as_str().unwrap()),
    )
    .unwrap();
    assert!(reopen(inside.path(), &name(&export))
        .unwrap_err()
        .to_string()
        .contains("outside"));
}
#[test]
fn reopen_count_limit_release_and_bounded_previews() {
    let (dir, _, export) = fixture("id\n001\n002\n", json!({}));
    let mut store = DatasetStore::new(dir.path(), dir.path()).unwrap();
    assert!(store
        .reopen(ReopenRequest {
            record_path: name(&export),
            preview_rows: 51
        })
        .is_err());
    let mut id = String::new();
    for _ in 0..16 {
        let opened = store
            .reopen(ReopenRequest {
                record_path: name(&export),
                preview_rows: 0,
            })
            .unwrap();
        assert_eq!(opened["row_count"], 2);
        assert_eq!(opened["preview"]["returned_rows"], 0);
        id = opened["dataset_id"].as_str().unwrap().into();
    }
    assert!(store
        .reopen(ReopenRequest {
            record_path: name(&export),
            preview_rows: 0
        })
        .is_err());
    store
        .release(serde_json::from_value(json!({"dataset_id":id})).unwrap())
        .unwrap();
    assert!(store
        .reopen(ReopenRequest {
            record_path: name(&export),
            preview_rows: 1
        })
        .is_ok());
    assert!(Path::new(export["record_path"].as_str().unwrap()).exists());
}
#[cfg(unix)]
#[test]
fn in_root_record_alias_resolves_ipc_beside_canonical_record() {
    use std::os::unix::fs::symlink;
    let (dir, _, export) = fixture("id\n001\n", json!({}));
    fs::create_dir(dir.path().join("aliases")).unwrap();
    symlink(
        export["record_path"].as_str().unwrap(),
        dir.path().join("aliases/record.json"),
    )
    .unwrap();
    let opened = reopen(dir.path(), "aliases/record.json").unwrap();
    assert_eq!(opened["preview"]["rows"], json!([["001"]]));
    assert_eq!(
        opened["reopen_verification"]["record_path"],
        "aliases/record.json"
    );
}
#[test]
fn reopened_payloads_charge_existing_retained_datasets_before_decode() {
    let long = "a".repeat(5 * 1024 * 1024);
    let (dir, _, export) = fixture(&format!("text\n{long}\n"), json!({}));
    let mut store = DatasetStore::new(dir.path(), dir.path()).unwrap();
    for _ in 0..12 {
        store
            .reopen(ReopenRequest {
                record_path: name(&export),
                preview_rows: 0,
            })
            .unwrap();
    }
    let error = store
        .reopen(ReopenRequest {
            record_path: name(&export),
            preview_rows: 0,
        })
        .unwrap_err()
        .to_string();
    assert!(error.contains("decoded allocation"), "{error}");
}

#[test]
fn escaped_wide_schema_and_bounded_provenance_export_records_over_32_kib_reopen() {
    let dir = tempfile::tempdir().unwrap();
    let names = (0..64)
        .map(|index| format!("c{index:02}{}", "\"".repeat(125)))
        .collect::<Vec<_>>();
    assert!(names.iter().all(|name| name.len() == 128));
    let mut csv = csv::Writer::from_writer(Vec::new());
    csv.write_record(&names).unwrap();
    csv.write_record(vec!["001"; 64]).unwrap();
    fs::write(dir.path().join("table.csv"), csv.into_inner().unwrap()).unwrap();
    let declared = "\"".repeat(256);
    let mut store = DatasetStore::new(dir.path(), dir.path()).unwrap();
    let mut current = store.open(serde_json::from_value(json!({
        "path":"table.csv", "preview_rows":0,
        "metadata":{"species":declared,"reference":declared,"coordinates":declared,"units":declared}
    })).unwrap()).unwrap();
    for _ in 0..31 {
        let previous = current["dataset_id"].as_str().unwrap().to_owned();
        current = store
            .query(serde_json::from_value(json!({"dataset_id":previous,"preview_rows":0})).unwrap())
            .unwrap();
        store
            .release(serde_json::from_value(json!({"dataset_id":previous})).unwrap())
            .unwrap();
    }
    let earlier = store
        .export(ExportRequest {
            dataset_id: current["dataset_id"].as_str().unwrap().into(),
        })
        .unwrap();
    let mut next_provenance = earlier["record"]["provenance"].clone();
    let filter = json!({"column":names[0],"op":"eq","value":""});
    next_provenance["operations"]
        .as_array_mut()
        .unwrap()
        .push(json!({
            "parent_dataset_id":current["dataset_id"], "filter":filter, "sort":null, "select":null,
            "execution_order":["filter","stable_sort_nulls_last","select"]
        }));
    let padding = 8190 - serde_json::to_vec(&next_provenance).unwrap().len();
    let final_dataset = store
        .query(
            serde_json::from_value(json!({
                "dataset_id":current["dataset_id"], "preview_rows":0,
                "filter":{"column":names[0],"op":"eq","value":"x".repeat(padding)}
            }))
            .unwrap(),
        )
        .unwrap();
    let exported = store
        .export(ExportRequest {
            dataset_id: final_dataset["dataset_id"].as_str().unwrap().into(),
        })
        .unwrap();
    assert_eq!(
        serde_json::to_vec(&exported["record"]["provenance"])
            .unwrap()
            .len(),
        8190
    );
    let record_length = fs::metadata(exported["record_path"].as_str().unwrap())
        .unwrap()
        .len();
    assert!(
        record_length > 32 * 1024,
        "test record was only {record_length} bytes"
    );
    assert!(record_length <= 64 * 1024);
    // Escaped 64-column names and near-limit provenance also fit both portable
    // companions without weakening the strict record's independent bound.
    for (field, limit) in [("manifest_path", 128 * 1024), ("readme_path", 16 * 1024)] {
        assert!(
            fs::metadata(exported[field].as_str().unwrap())
                .unwrap()
                .len()
                <= limit
        );
    }
    let recovered = reopen(dir.path(), &name(&exported)).unwrap();
    assert_eq!(recovered["schema"], final_dataset["schema"]);
    assert_eq!(recovered["row_count"], final_dataset["row_count"]);
    assert_eq!(
        recovered["scientific_metadata"],
        final_dataset["scientific_metadata"]
    );
}
