//! One established CSV decoder for both allocation preflight and materialization.
//! `csv-core` is the engine used by `csv`; Polars never reparses the source.
use super::{checked, err, memory_charge, ColumnType, Result, MAX_COLUMNS};
use csv_core::{ReadFieldResult, Reader};
use polars::prelude::*;
use std::collections::{BTreeMap, HashSet};

struct Fields<'a> {
    reader: Reader,
    remaining: &'a [u8],
    first_field: bool,
    // A decoded field cannot exceed the bounded source snapshot. This scratch
    // buffer is temporary parser memory, not part of the retained-table charge.
    output: Vec<u8>,
}
struct Field<'a> {
    value: &'a str,
    quoted_empty: bool,
    record_end: bool,
}
impl<'a> Fields<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self {
            reader: Reader::new(),
            remaining: bytes,
            first_field: true,
            output: vec![0; bytes.len().max(1)],
        }
    }
    fn next(&mut self) -> Result<Option<Field<'_>>> {
        let mut written = 0;
        let raw_start = self.remaining;
        loop {
            let (result, consumed, produced) = self
                .reader
                .read_field(self.remaining, &mut self.output[written..]);
            self.remaining = &self.remaining[consumed..];
            written += produced;
            match result {
                ReadFieldResult::Field { record_end } => {
                    let value = checked(std::str::from_utf8(&self.output[..written]))?;
                    let raw = &raw_start[..raw_start.len() - self.remaining.len()];
                    let quoted = validate_quote_envelope(raw, value, record_end, self.first_field)?;
                    self.first_field = false;
                    return Ok(Some(Field {
                        value,
                        quoted_empty: written == 0 && quoted,
                        record_end,
                    }));
                }
                ReadFieldResult::End => return Ok(None),
                ReadFieldResult::InputEmpty => continue, // Flush the final field at EOF.
                ReadFieldResult::OutputFull => {
                    return Err(err("CSV field exceeds source byte bound"))
                }
            }
        }
    }
    fn headers(&mut self) -> Result<Vec<String>> {
        let mut headers = Vec::new();
        while let Some(field) = self.next()? {
            if headers.len() == MAX_COLUMNS {
                return Err(err("table must have 1..64 columns"));
            }
            if field.value.is_empty()
                || field.value.len() > 128
                || field.value.chars().any(char::is_control)
            {
                return Err(err(
                    "column names must be 1..128 bytes without control characters",
                ));
            }
            headers.push(field.value.to_owned());
            if field.record_end {
                break;
            }
        }
        if headers.is_empty() {
            return Err(err("table must have 1..64 columns"));
        }
        if headers.iter().collect::<HashSet<_>>().len() != headers.len() {
            return Err(err("duplicate CSV column names are not permitted"));
        }
        Ok(headers)
    }
    fn records(
        &mut self,
        width: usize,
        mut visit: impl FnMut(usize, usize, Field<'_>) -> Result<()>,
    ) -> Result<usize> {
        let (mut rows, mut column) = (0, 0);
        while let Some(field) = self.next()? {
            if column >= width {
                return Err(err("CSV row width does not match header"));
            }
            let record_end = field.record_end;
            visit(rows, column, field)?;
            column += 1;
            if record_end {
                if column != width {
                    return Err(err("CSV row width does not match header"));
                }
                rows += 1;
                column = 0;
            }
        }
        if column != 0 {
            return Err(err("incomplete CSV record"));
        }
        Ok(rows)
    }
}

// csv-core determines all field/record boundaries and unescapes content, but
// deliberately accepts broken quoting. Validate its consumed field's envelope
// using the standard doubled-quote length identity, without a second CSV parser.
fn validate_quote_envelope(raw: &[u8], value: &str, record_end: bool, first: bool) -> Result<bool> {
    let raw = if first {
        raw.strip_prefix(b"\xef\xbb\xbf").unwrap_or(raw)
    } else {
        raw
    };
    // csv-core can consume ignored physical blank lines before the next field.
    let raw = &raw[raw
        .iter()
        .position(|byte| !matches!(byte, b'\r' | b'\n'))
        .unwrap_or(raw.len())..];
    let raw = if record_end {
        raw.strip_suffix(b"\r")
            .or_else(|| raw.strip_suffix(b"\n"))
            .unwrap_or(raw)
    } else {
        raw.strip_suffix(b",").unwrap_or(raw)
    };
    let quoted = raw.starts_with(b"\"");
    if quoted
        && (!raw.ends_with(b"\"")
            || raw.len() != value.len() + value.bytes().filter(|byte| *byte == b'"').count() + 2)
    {
        return Err(err(
            "invalid CSV quoting: unterminated field or characters after closing quote",
        ));
    }
    Ok(quoted)
}

enum Cell<'a> {
    String(Option<&'a str>),
    Int64(Option<i64>),
    Float64(Option<f64>),
    Boolean(Option<bool>),
}
fn cell<'a>(kind: &ColumnType, field: Field<'a>) -> Result<Cell<'a>> {
    let value = field.value;
    Ok(match kind {
        ColumnType::String => Cell::String(if value.is_empty() && !field.quoted_empty {
            None
        } else {
            Some(value)
        }),
        ColumnType::Int64 => {
            let value = value.trim_start_matches([' ', '\t']);
            Cell::Int64(if value.is_empty() {
                None
            } else {
                Some(checked(value.parse::<i64>())?)
            })
        }
        ColumnType::Float64 => {
            let value = value.trim_start_matches([' ', '\t']);
            let number = if value.is_empty() {
                None
            } else {
                Some(checked(value.parse::<f64>())?)
            };
            if number.is_some_and(|n| !n.is_finite()) {
                return Err(err("non-finite float64 values are not supported; preserve them as strings or clean explicitly"));
            }
            Cell::Float64(number)
        }
        ColumnType::Boolean => Cell::Boolean(if value.is_empty() {
            None
        } else if value.eq_ignore_ascii_case("true") {
            Some(true)
        } else if value.eq_ignore_ascii_case("false") {
            Some(false)
        } else {
            return Err(err("invalid boolean: expected true or false"));
        }),
    })
}

struct Plan {
    headers: Vec<String>,
    kinds: Vec<ColumnType>,
    rows: usize,
    charge: usize,
}
impl Plan {
    // No Polars column, row vector or table builder is allocated in preflight.
    fn inspect(
        bytes: &[u8],
        declared: &BTreeMap<String, ColumnType>,
        available: usize,
    ) -> Result<Self> {
        let mut fields = Fields::new(bytes);
        let headers = fields.headers()?;
        if declared.keys().any(|name| !headers.contains(name)) {
            return Err(err("schema names must match existing CSV columns"));
        }
        let kinds: Vec<_> = headers
            .iter()
            .map(|name| declared.get(name).cloned().unwrap_or(ColumnType::String))
            .collect();
        // One validity bitmap per column and an additional value bitmap for
        // Boolean. Charging even absent validity buffers is conservative.
        let bitmaps = kinds.len()
            + kinds
                .iter()
                .filter(|kind| matches!(kind, ColumnType::Boolean))
                .count();
        let (mut cells_and_payload, mut charge) = (0usize, 0usize);
        let rows = fields.records(headers.len(), |row, column, field| {
            let payload = match cell(&kinds[column], field)? {
                Cell::String(value) => value.map_or(0, str::len),
                Cell::Int64(_) | Cell::Float64(_) => 8,
                Cell::Boolean(_) => 0,
            };
            cells_and_payload = cells_and_payload.saturating_add(17).saturating_add(payload);
            charge =
                cells_and_payload.saturating_add((row + 1).div_ceil(8).saturating_mul(bitmaps));
            if charge > available {
                return Err(err(
                    "input row/column allocation exceeds remaining retained dataset budget",
                ));
            }
            Ok(())
        })?;
        Ok(Self {
            headers,
            kinds,
            rows,
            charge,
        })
    }
}

enum Builder {
    String(StringChunkedBuilder),
    Int64(PrimitiveChunkedBuilder<Int64Type>),
    Float64(PrimitiveChunkedBuilder<Float64Type>),
    Boolean(BooleanChunkedBuilder),
}
impl Builder {
    fn new(name: &str, kind: &ColumnType, rows: usize) -> Self {
        match kind {
            ColumnType::String => Self::String(StringChunkedBuilder::new(name.into(), rows)),
            ColumnType::Int64 => Self::Int64(PrimitiveChunkedBuilder::new(name.into(), rows)),
            ColumnType::Float64 => Self::Float64(PrimitiveChunkedBuilder::new(name.into(), rows)),
            ColumnType::Boolean => Self::Boolean(BooleanChunkedBuilder::new(name.into(), rows)),
        }
    }
    fn append(&mut self, value: Cell<'_>) -> Result<()> {
        match (self, value) {
            (Self::String(builder), Cell::String(value)) => builder.append_option(value),
            (Self::Int64(builder), Cell::Int64(value)) => builder.append_option(value),
            (Self::Float64(builder), Cell::Float64(value)) => builder.append_option(value),
            (Self::Boolean(builder), Cell::Boolean(value)) => builder.append_option(value),
            _ => return Err(err("CSV type changed after preflight")),
        }
        Ok(())
    }
    fn finish(self) -> Column {
        match self {
            Self::String(builder) => builder.finish().into_series().into(),
            Self::Int64(builder) => builder.finish().into_series().into(),
            Self::Float64(builder) => builder.finish().into_series().into(),
            Self::Boolean(builder) => builder.finish().into_series().into(),
        }
    }
}

pub(super) fn read(
    bytes: &[u8],
    declared: &BTreeMap<String, ColumnType>,
    available: usize,
) -> Result<DataFrame> {
    let plan = Plan::inspect(bytes, declared, available)?;
    // Only allocate typed columns after the entire snapshot has passed the
    // authoritative parser's row, type and remaining-retained-budget checks.
    let mut builders: Vec<_> = plan
        .headers
        .iter()
        .zip(&plan.kinds)
        .map(|(name, kind)| Builder::new(name, kind, plan.rows))
        .collect();
    let mut fields = Fields::new(bytes);
    if fields.headers()? != plan.headers {
        return Err(err("CSV header changed after preflight"));
    }
    let rows = fields.records(plan.headers.len(), |_, column, field| {
        builders[column].append(cell(&plan.kinds[column], field)?)
    })?;
    let frame = checked(DataFrame::new(
        builders.into_iter().map(Builder::finish).collect(),
    ))?;
    // Defense in depth, not the allocation guard: both passes used the same
    // decoder over identical bytes and builders received exactly its records.
    if rows != plan.rows || frame.height() != plan.rows || memory_charge(&frame) > plan.charge {
        return Err(err(
            "CSV materialization disagrees with validated rows or allocation bound",
        ));
    }
    Ok(frame)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preflight_includes_payload_bitmaps_and_exact_remaining_budget() {
        let bytes = b"text,count,ratio,active\r\nlong decoded payload,9007199254740993,1.5,true\r\n\"\",,,\r\n";
        let declared = BTreeMap::from([
            ("count".into(), ColumnType::Int64),
            ("ratio".into(), ColumnType::Float64),
            ("active".into(), ColumnType::Boolean),
        ]);
        let plan = Plan::inspect(bytes, &declared, usize::MAX).unwrap();
        assert_eq!(plan.rows, 2);
        // Eight cells, two 8-byte columns, string payload, five packed bitmaps.
        assert_eq!(
            plan.charge,
            8 * 17 + 2 * 2 * 8 + "long decoded payload".len() + 5
        );
        assert!(Plan::inspect(bytes, &declared, plan.charge).is_ok());
        assert!(Plan::inspect(bytes, &declared, plan.charge - 1).is_err());
        assert!(read(bytes, &declared, plan.charge - 1).is_err());
        let frame = read(bytes, &declared, plan.charge).unwrap();
        assert!(memory_charge(&frame) <= plan.charge);
        // The caller passes 64 MiB minus all already-retained tables.
        let retained = super::super::MAX_MEMORY_BYTES - plan.charge + 1;
        assert!(
            Plan::inspect(bytes, &declared, super::super::MAX_MEMORY_BYTES - retained).is_err()
        );
    }

    #[test]
    fn cell_budget_rejects_typed_payload_before_any_column_builder() {
        let names: Vec<_> = (0..64).map(|i| format!("c{i}")).collect();
        let declared = names
            .iter()
            .map(|name| (name.clone(), ColumnType::Int64))
            .collect();
        // The former 17-byte cell-only guard accepted this ~5 MiB input,
        // but the complete int64 payload makes it exceed 64 MiB.
        let csv = format!(
            "{}\n{}",
            names.join(","),
            format!("{}\n", vec!["1"; 64].join(",")).repeat(42_000)
        );
        assert!(csv.len() < super::super::MAX_INPUT_BYTES as usize);
        assert!(42_000 * names.len() * 17 < super::super::MAX_MEMORY_BYTES);
        let error = Plan::inspect(csv.as_bytes(), &declared, super::super::MAX_MEMORY_BYTES)
            .err()
            .unwrap();
        assert!(error.to_string().contains("allocation"));
    }

    #[test]
    fn blank_lines_never_become_allocated_rows() {
        let names: Vec<_> = (0..64).map(|i| format!("c{i}")).collect();
        let csv = format!("{}\r\n{}", names.join(","), "\r\n\n\r".repeat(62_000));
        let plan = Plan::inspect(csv.as_bytes(), &BTreeMap::new(), 0).unwrap();
        assert_eq!((plan.rows, plan.charge), (0, 0));
        assert_eq!(
            read(csv.as_bytes(), &BTreeMap::new(), 0).unwrap().shape(),
            (0, 64)
        );
    }

    #[test]
    fn authoritative_decoder_agrees_with_csv_records_for_generated_fields() {
        let values = [
            "",
            "001",
            "comma,field",
            "a\"b",
            "\n",
            "\r",
            "\r\n",
            "a\n\nb",
            "λ🧬",
        ];
        for terminator in [
            csv::Terminator::CRLF,
            csv::Terminator::Any(b'\n'),
            csv::Terminator::Any(b'\r'),
        ] {
            let mut writer = csv::WriterBuilder::new()
                .terminator(terminator)
                .from_writer(Vec::new());
            writer.write_record(["a", "b"]).unwrap();
            for left in &values {
                for right in &values {
                    writer.write_record([left, right]).unwrap();
                }
            }
            let bytes = writer.into_inner().unwrap();
            let mut expected = csv::Reader::from_reader(bytes.as_slice());
            let expected: Vec<_> = expected
                .records()
                .map(std::result::Result::unwrap)
                .collect();
            let mut fields = Fields::new(&bytes);
            assert_eq!(fields.headers().unwrap(), ["a", "b"]);
            let rows = fields
                .records(2, |row, column, field| {
                    assert_eq!(field.value, &expected[row][column]);
                    Ok(())
                })
                .unwrap();
            assert_eq!(rows, values.len() * values.len());
        }
    }
}
