//! Regression checks added during independent contract review.

use biov_data::{DatasetStore, OpenRequest, QueryRequest};
use serde_json::{json, Value};
use std::fs;

fn fixture(csv: &str) -> (tempfile::TempDir, DatasetStore) {
    let directory = tempfile::tempdir().unwrap();
    let inputs = directory.path().join("inputs");
    let outputs = directory.path().join("outputs");
    fs::create_dir(&inputs).unwrap();
    fs::create_dir(&outputs).unwrap();
    fs::write(inputs.join("table.csv"), csv).unwrap();
    let store = DatasetStore::new(&inputs, &outputs).unwrap();
    (directory, store)
}

fn query(store: &mut DatasetStore, value: Value) -> Value {
    let request: QueryRequest = serde_json::from_value(value).unwrap();
    store.query(request).unwrap()
}

#[test]
fn boolean_equality_keeps_true_false_and_null_distinct() {
    let (_directory, mut store) = fixture(
        "label,passed\nfirst_true,true\nfirst_false,false\nmissing,\nsecond_true,true\nsecond_false,false\n",
    );
    let request: OpenRequest = serde_json::from_value(json!({
        "path": "table.csv",
        "schema": {"passed": "boolean"},
        "preview_rows": 1
    }))
    .unwrap();
    let source = store.open(request).unwrap();
    assert_eq!(source["row_count"], 5);
    let id = source["dataset_id"].as_str().unwrap();

    let true_rows = query(
        &mut store,
        json!({
            "dataset_id": id,
            "filter": {"column": "passed", "op": "eq", "value": true}
        }),
    );
    assert_eq!(true_rows["row_count"], 2);
    assert_eq!(
        true_rows["preview"]["rows"],
        json!([["first_true", true], ["second_true", true]])
    );

    let false_rows = query(
        &mut store,
        json!({
            "dataset_id": id,
            "filter": {"column": "passed", "op": "eq", "value": false}
        }),
    );
    assert_eq!(false_rows["row_count"], 2);
    assert_eq!(
        false_rows["preview"]["rows"],
        json!([["first_false", false], ["second_false", false]])
    );

    let null_rows = query(
        &mut store,
        json!({
            "dataset_id": id,
            "filter": {"column": "passed", "op": "is_null"}
        }),
    );
    assert_eq!(null_rows["row_count"], 1);
    assert_eq!(null_rows["preview"]["rows"], json!([["missing", null]]));
}

#[test]
fn boolean_filter_requires_boolean_values_and_supported_operators() {
    let (_directory, mut store) = fixture("label,passed\na,true\nb,false\nc,\n");
    let source = store
        .open(
            serde_json::from_value(json!({
                "path": "table.csv", "schema": {"passed": "boolean"}
            }))
            .unwrap(),
        )
        .unwrap();
    let id = source["dataset_id"].as_str().unwrap();
    for filter in [
        json!({"column": "passed", "op": "eq", "value": "true"}),
        json!({"column": "passed", "op": "eq", "value": 1}),
        json!({"column": "passed", "op": "eq", "value": null}),
        json!({"column": "passed", "op": "gt", "value": true}),
        json!({"column": "passed", "op": "is_null", "value": false}),
    ] {
        let request = serde_json::from_value(json!({
            "dataset_id": id, "filter": filter
        }))
        .unwrap();
        assert!(store.query(request).is_err());
    }
}

#[cfg(unix)]
#[test]
fn non_utf8_roots_are_rejected_before_any_dataset_or_artifact_operation() {
    use std::{ffi::OsString, os::unix::ffi::OsStringExt};

    let directory = tempfile::tempdir().unwrap();
    let normal = directory.path().join("normal");
    let non_utf8 = directory
        .path()
        .join(OsString::from_vec(vec![b'd', b'i', b'r', 0xff]));
    fs::create_dir(&normal).unwrap();
    fs::create_dir(&non_utf8).unwrap();

    for (input, output) in [(&non_utf8, &normal), (&normal, &non_utf8)] {
        let result = DatasetStore::new(input, output);
        assert!(result.is_err());
        assert!(result.err().unwrap().to_string().contains("UTF-8"));
    }
    assert_eq!(fs::read_dir(&normal).unwrap().count(), 0);
    assert_eq!(fs::read_dir(&non_utf8).unwrap().count(), 0);
}

#[test]
fn artifact_and_derivation_limits_fail_without_losing_prior_outputs() {
    let (directory, mut store) = fixture("label\na\n");
    let mut current = store
        .open(serde_json::from_value(json!({"path":"table.csv"})).unwrap())
        .unwrap();
    for _ in 0..32 {
        let id = current["dataset_id"].as_str().unwrap().to_owned();
        current = query(&mut store, json!({"dataset_id":id}));
        store
            .release(serde_json::from_value(json!({"dataset_id":id})).unwrap())
            .unwrap();
    }
    assert!(store
        .query(serde_json::from_value(json!({"dataset_id":current["dataset_id"]})).unwrap())
        .is_err());
    for _ in 0..64 {
        store
            .export(serde_json::from_value(json!({"dataset_id":current["dataset_id"]})).unwrap())
            .unwrap();
    }
    assert!(store
        .export(serde_json::from_value(json!({"dataset_id":current["dataset_id"]})).unwrap())
        .is_err());
    assert_eq!(
        fs::read_dir(directory.path().join("outputs"))
            .unwrap()
            .count(),
        256 // 64 portable bundles, each with Arrow, record, manifest and README.
    );
}

#[test]
fn quoted_empty_string_remains_distinct_from_missing_field() {
    let (_directory, mut store) = fixture("label,value\nquoted,\"\"\nmissing,\n");
    let opened = store
        .open(serde_json::from_value(json!({"path":"table.csv"})).unwrap())
        .unwrap();
    assert_eq!(
        opened["preview"]["rows"],
        json!([["quoted", ""], ["missing", null]])
    );
}

#[test]
fn excessive_cell_allocation_is_rejected_before_polars_materialization() {
    let headers = (0..64)
        .map(|i| format!("c{i}"))
        .collect::<Vec<_>>()
        .join(",");
    let row = format!("{}\n", vec!["x"; 64].join(","));
    let (_directory, mut store) = fixture(&format!("{headers}\n{}", row.repeat(62000)));
    let result = store.open(serde_json::from_value(json!({"path":"table.csv"})).unwrap());
    assert!(result.unwrap_err().to_string().contains("allocation"));
}
