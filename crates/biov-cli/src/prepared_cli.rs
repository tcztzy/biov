//! Thin prepared-artifact adapter, sharing the bounded one-shot JSON transport.
use std::path::Path;

use biov_prepared::{PrepareFastaRequest, PreparedStore};

use crate::storage_cli::{read_request, write_result};

pub(crate) fn fasta(store_root: &Path, request_file: &Path) -> Result<(), String> {
    let request: PrepareFastaRequest = read_request(request_file)?;
    let store = PreparedStore::new(store_root).map_err(|error| error.to_string())?;
    let result = store
        .prepare_fasta(request)
        .map_err(|error| error.to_string())?;
    write_result(&result)
}
