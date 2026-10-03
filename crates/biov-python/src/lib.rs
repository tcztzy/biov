//! Transport and exception adaptation only; biological rules live in biov-core.
use biov_core::sequence::{self, Kind, SequenceError};
use pyo3::{
    create_exception,
    exceptions::{PyTypeError, PyUnicodeEncodeError, PyValueError},
    prelude::*,
    types::PyList,
};

create_exception!(biov._native, SequenceValidationError, PyValueError);

fn map_error(error: SequenceError) -> PyErr {
    match error {
        SequenceError::InvalidSymbols { .. } => SequenceValidationError::new_err(error.to_string()),
        SequenceError::NotNucleic => PyTypeError::new_err(error.to_string()),
        SequenceError::UnknownKind(_) => PyValueError::new_err(error.to_string()),
    }
}

fn extract_sequences(py: Python<'_>, values: &Bound<'_, PyList>) -> PyResult<Vec<Option<String>>> {
    values.extract().map_err(|error: PyErr| {
        if error.is_instance_of::<PyUnicodeEncodeError>(py) {
            SequenceValidationError::new_err("invalid sequence symbol: unpaired Unicode surrogate")
        } else {
            error
        }
    })
}

#[pyfunction]
#[pyo3(signature = (values, *, kind))]
fn normalize_sequences(
    py: Python<'_>,
    values: &Bound<'_, PyList>,
    kind: &str,
) -> PyResult<Vec<Option<String>>> {
    let kind: Kind = kind.parse().map_err(map_error)?;
    let values = extract_sequences(py, values)?;
    py.detach(move || sequence::normalize_batch(&values, kind))
        .map_err(map_error)
}

#[pyfunction]
#[pyo3(signature = (values, *, kind))]
fn reverse_complements(
    py: Python<'_>,
    values: &Bound<'_, PyList>,
    kind: &str,
) -> PyResult<Vec<Option<String>>> {
    let kind: Kind = kind.parse().map_err(map_error)?;
    let values = extract_sequences(py, values)?;
    py.detach(move || sequence::reverse_complement_batch(&values, kind))
        .map_err(map_error)
}

/// Return symbol counts for nullable DNA, RNA or protein strings after validation.
#[pyfunction]
#[pyo3(signature = (values, *, kind))]
fn sequence_lengths(
    py: Python<'_>,
    values: &Bound<'_, PyList>,
    kind: &str,
) -> PyResult<Vec<Option<usize>>> {
    let kind: Kind = kind.parse().map_err(map_error)?;
    let values = extract_sequences(py, values)?;
    py.detach(move || sequence::length_batch(&values, kind))
        .map_err(map_error)
}

/// Return mean IUPAC GC probabilities for DNA or RNA; empty sequences yield 0.0.
#[pyfunction]
#[pyo3(signature = (values, *, kind))]
fn weighted_gc_fractions(
    py: Python<'_>,
    values: &Bound<'_, PyList>,
    kind: &str,
) -> PyResult<Vec<Option<f64>>> {
    let kind: Kind = kind.parse().map_err(map_error)?;
    let values = extract_sequences(py, values)?;
    py.detach(move || sequence::weighted_gc_fraction_batch(&values, kind))
        .map_err(map_error)
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "SequenceValidationError",
        m.py().get_type::<SequenceValidationError>(),
    )?;
    m.add_function(wrap_pyfunction!(normalize_sequences, m)?)?;
    m.add_function(wrap_pyfunction!(reverse_complements, m)?)?;
    m.add_function(wrap_pyfunction!(sequence_lengths, m)?)?;
    m.add_function(wrap_pyfunction!(weighted_gc_fractions, m)?)?;
    Ok(())
}
