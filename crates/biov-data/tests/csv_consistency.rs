use biov_data::{DatasetStore, ExportRequest};
use polars::prelude::*;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{fs, io::Cursor};

fn check(csv: &str, expected: Value) {
    let directory = tempfile::tempdir().unwrap();
    fs::write(directory.path().join("table.csv"), csv).unwrap();
    let mut store = DatasetStore::new(directory.path(), directory.path()).unwrap();
    let opened = store
        .open(serde_json::from_value(json!({"path":"table.csv", "preview_rows":50})).unwrap())
        .unwrap();
    assert_eq!(
        opened["row_count"],
        expected.as_array().unwrap().len(),
        "{csv:?}"
    );
    assert_eq!(opened["preview"]["rows"], expected, "{csv:?}");
    assert_eq!(
        opened["source_sha256"],
        format!("{:x}", Sha256::digest(csv.as_bytes()))
    );
    let exported = store
        .export(ExportRequest {
            dataset_id: opened["dataset_id"].as_str().unwrap().into(),
        })
        .unwrap();
    assert_eq!(exported["record"]["row_count"], opened["row_count"]);
    let bytes = fs::read(exported["execution_host_path"].as_str().unwrap()).unwrap();
    let frame = IpcReader::new(Cursor::new(bytes)).finish().unwrap();
    assert_eq!(frame.height(), expected.as_array().unwrap().len());
    for (row_index, row) in expected.as_array().unwrap().iter().enumerate() {
        for (column, value) in row.as_array().unwrap().iter().enumerate() {
            assert_eq!(
                frame.get_columns()[column].str().unwrap().get(row_index),
                value.as_str()
            );
        }
    }
}

#[test]
fn lf_crlf_cr_and_mixed_terminators_preserve_complete_records() {
    for terminator in ["\n", "\r\n", "\r"] {
        check(
            &format!("x{terminator}a{terminator}b{terminator}"),
            json!([["a"], ["b"]]),
        );
        check(
            &format!(
                "{terminator}x,y{terminator}a,b{terminator}{terminator}c,d{terminator}{terminator}"
            ),
            json!([["a", "b"], ["c", "d"]]),
        );
    }
    check("\r\nx,y\ra,b\n\r\nc,d\r\n", json!([["a", "b"], ["c", "d"]]));
}

#[test]
fn quoted_embedded_terminators_blank_lines_escaped_quotes_and_nulls_survive() {
    check(
        "x,y\n\"a\n\nb\",\"c\rd\"\r\n\"a\r\n\r\nb\",\"say \"\"hello\"\"\"\r,\"\"\n",
        json!([
            ["a\n\nb", "c\rd"],
            ["a\r\n\r\nb", "say \"hello\""],
            [null, ""]
        ]),
    );
    check("x\n\"\"\n\n", json!([[""]]));
    check("x,y\na,\nb,\"\"\n", json!([["a", null], ["b", ""]]));
    check("x,y\n,\n\"\",\"\"", json!([[null, null], ["", ""]]));
    check(
        "\u{feff}x,y\r001,λ🧬\r002,",
        json!([["001", "λ🧬"], ["002", null]]),
    );
}

#[test]
fn header_only_and_missing_final_record_terminator_are_complete() {
    check("x,y", json!([]));
    check("x,y\r\n\r\n", json!([]));
    check("x,y\ra,b", json!([["a", "b"]]));
    check("x,y\na,\"\"", json!([["a", ""]]));
}

#[test]
fn ragged_rows_are_rejected_for_every_record_terminator() {
    let directory = tempfile::tempdir().unwrap();
    let mut store = DatasetStore::new(directory.path(), directory.path()).unwrap();
    for terminator in ["\n", "\r\n", "\r"] {
        for row in ["a", "a,b,c", "\"\""] {
            fs::write(
                directory.path().join("table.csv"),
                format!("x,y{terminator}{row}{terminator}"),
            )
            .unwrap();
            let result = store.open(serde_json::from_value(json!({"path":"table.csv"})).unwrap());
            assert!(result.unwrap_err().to_string().contains("width"));
        }
    }
}

#[test]
fn nullable_explicit_types_preserve_int64_precision_and_string_spelling() {
    let directory = tempfile::tempdir().unwrap();
    fs::write(directory.path().join("table.csv"), "id,count,ratio,active\r001,9007199254740993,1.25,TRUE\r\n\r002,-9223372036854775808,-0.0,false\n003,,,\r004,\"\",\"\",\"\"\r").unwrap();
    let mut store = DatasetStore::new(directory.path(), directory.path()).unwrap();
    let opened = store.open(serde_json::from_value(json!({"path":"table.csv", "schema":{"count":"int64", "ratio":"float64", "active":"boolean"}})).unwrap()).unwrap();
    assert_eq!(opened["row_count"], 4);
    assert_eq!(
        opened["preview"]["rows"],
        json!([
            ["001", 9007199254740993i64, 1.25, true],
            ["002", i64::MIN, -0.0, false],
            ["003", null, null, null],
            ["004", null, null, null]
        ])
    );
}

#[test]
fn damaged_quoted_fields_fail_instead_of_silently_changing_content() {
    let directory = tempfile::tempdir().unwrap();
    let mut store = DatasetStore::new(directory.path(), directory.path()).unwrap();
    for csv in [
        "x\n\"a\n",
        "x\n\"a\"junk\n",
        "x\n\"a\"junk\"\n",
        "x\n\"",
        "x\n\"\"\"",
        "x,y\na,\"broken\r\n",
        "\"broken\nheader\n",
    ] {
        fs::write(directory.path().join("table.csv"), csv).unwrap();
        let error = store
            .open(serde_json::from_value(json!({"path":"table.csv"})).unwrap())
            .unwrap_err();
        assert!(error.to_string().contains("quoting"), "{csv:?}: {error}");
    }
}

#[test]
fn numeric_empty_and_leading_whitespace_policy_is_preserved() {
    let directory = tempfile::tempdir().unwrap();
    let mut store = DatasetStore::new(directory.path(), directory.path()).unwrap();
    fs::write(
        directory.path().join("table.csv"),
        "x\n\"\"\n \n\" \"\n\t\n+1\n 1\n\" 1\"\n",
    )
    .unwrap();
    let opened = store
        .open(
            serde_json::from_value(
                json!({"path":"table.csv", "schema":{"x":"int64"},"preview_rows":50}),
            )
            .unwrap(),
        )
        .unwrap();
    assert_eq!(
        opened["preview"]["rows"],
        json!([[null], [null], [null], [null], [1], [1], [1]])
    );
    fs::write(directory.path().join("table.csv"), "x\n1 \n").unwrap();
    assert!(store
        .open(serde_json::from_value(json!({"path":"table.csv", "schema":{"x":"int64"}})).unwrap())
        .is_err());
}
