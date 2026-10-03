//! Bounded input policy around upstream FASTA parsing and FAI generation.
//! This module does not decode identifiers, calculate offsets or validate wrapping.
use crate::*;
use noodles_fasta::{fai, io::Indexer};
use sha2::{Digest, Sha256};
use std::{
    collections::HashSet,
    io::{self, BufRead, BufReader, Read, Write},
};

/// Presents complete, bounded physical lines, avoiding upstream chunk-boundary
/// ambiguity while bounding its definition allocation. Biological records and
/// index offsets are exclusively decoded by noodles-fasta.
struct BoundedLines<R> {
    inner: R,
    line: Vec<u8>,
    position: usize,
    digest: Sha256,
    bytes: u64,
}
impl<R: BufRead> BoundedLines<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            line: Vec::with_capacity(MAX_LINE_BYTES),
            position: 0,
            digest: Sha256::new(),
            bytes: 0,
        }
    }
    fn next_line(&mut self) -> io::Result<()> {
        self.line.clear();
        self.position = 0;
        loop {
            let chunk = self.inner.fill_buf()?;
            if chunk.is_empty() {
                break;
            }
            let n = chunk
                .iter()
                .position(|&b| b == b'\n')
                .map_or(chunk.len(), |i| i + 1);
            if n > MAX_LINE_BYTES - self.line.len() {
                return Err(io::Error::new(
                    io::ErrorKind::FileTooLarge,
                    "physical FASTA line exceeds 1 MiB",
                ));
            }
            let done = chunk[n - 1] == b'\n';
            self.line.extend_from_slice(&chunk[..n]);
            self.digest.update(&chunk[..n]);
            self.bytes = self
                .bytes
                .checked_add(n as u64)
                .ok_or_else(|| io::Error::other("source byte count overflow"))?;
            self.inner.consume(n);
            if done {
                break;
            }
        }
        if self.line.is_empty() {
            return Ok(());
        }
        let body = if let Some(body) = self.line.strip_suffix(b"\n") {
            body.strip_suffix(b"\r").unwrap_or(body)
        } else {
            &self.line
        };
        if body.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "empty physical lines are unsupported",
            ));
        }
        if body[0] == b'>' {
            if !body
                .iter()
                .all(|&b| (b' '..=b'~').contains(&b) || b == b'\t')
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "definitions must be printable ASCII (horizontal tab permitted)",
                ));
            }
        } else if !body
            .iter()
            .all(|b| b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(b))
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "sequence lines must contain only case-preserved IUPAC DNA letters",
            ));
        }
        Ok(())
    }
}
impl<R: BufRead> BufRead for BoundedLines<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.position == self.line.len() {
            self.next_line()?;
        }
        Ok(&self.line[self.position..])
    }
    fn consume(&mut self, amount: usize) {
        self.position += amount.min(self.line.len() - self.position);
    }
}
impl<R: BufRead> Read for BoundedLines<R> {
    fn read(&mut self, destination: &mut [u8]) -> io::Result<usize> {
        if destination.is_empty() {
            return Ok(0);
        }
        let data = self.fill_buf()?;
        let n = data.len().min(destination.len());
        destination[..n].copy_from_slice(&data[..n]);
        self.consume(n);
        Ok(n)
    }
}

struct HashWriter<W> {
    inner: W,
    digest: Sha256,
    bytes: u64,
}
impl<W: Write> HashWriter<W> {
    fn new(inner: W) -> Self {
        Self {
            inner,
            digest: Sha256::new(),
            bytes: 0,
        }
    }
    fn identity(self, path: &str) -> OutputIdentity {
        OutputIdentity {
            path: path.into(),
            bytes: self.bytes,
            sha256: format!("{:x}", self.digest.finalize()),
        }
    }
}
impl<W: Write> Write for HashWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if bytes.len() as u64 > MAX_METADATA_BYTES as u64 - self.bytes {
            return Err(io::Error::new(
                io::ErrorKind::FileTooLarge,
                "generated index or dictionary exceeds 16 MiB",
            ));
        }
        let n = self.inner.write(bytes)?;
        self.digest.update(&bytes[..n]);
        self.bytes += n as u64;
        Ok(n)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) struct Indexed {
    pub sequence_count: usize,
    pub total_bases: u64,
    pub source_sha256: String,
    pub source_bytes: u64,
    pub fai: OutputIdentity,
    pub dictionary: OutputIdentity,
}

pub(crate) fn index<R: Read, F: Write, D: Write>(
    source: R,
    fai_output: F,
    dictionary_output: D,
) -> Result<Indexed, PreparedError> {
    index_with_capacity(source, fai_output, dictionary_output, IO_BUFFER_BYTES)
}
fn index_with_capacity<R: Read, F: Write, D: Write>(
    source: R,
    fai_output: F,
    dictionary_output: D,
    capacity: usize,
) -> Result<Indexed, PreparedError> {
    let mut source = BoundedLines::new(BufReader::with_capacity(capacity, source));
    let mut indexer = Indexer::new(&mut source);
    let mut fai_output = HashWriter::new(fai_output);
    let mut dictionary_output = HashWriter::new(dictionary_output);
    dictionary_output
        .write_all(b"sequence_id\tlength\n")
        .map_err(index_io)?;
    let mut seen = HashSet::new();
    let mut seen_charge = 0usize;
    let mut sequence_count = 0usize;
    let mut total_bases = 0u64;
    while let Some(record) = indexer.index_record().map_err(|e| index_io(e.into()))? {
        if sequence_count == MAX_RECORDS {
            return Err(limit("FASTA contains more than 100000 records"));
        }
        let name: &[u8] = record.name().as_ref();
        if name.len() > MAX_ID_BYTES {
            return Err(limit("FASTA identifier exceeds 4096 bytes"));
        }
        if name.is_empty() || !name.iter().all(|&b| (b'!'..=b'~').contains(&b)) {
            return Err(PreparedError::InvalidFasta(
                "identifiers must be nonempty printable ASCII tokens".into(),
            ));
        }
        seen_charge = seen_charge
            .checked_add(name.len() + 128)
            .ok_or_else(|| limit("identifier metadata overflow"))?;
        if seen_charge > MAX_METADATA_BYTES {
            return Err(limit("identifier metadata charge exceeds 16 MiB"));
        }
        if !seen.insert(name.to_vec()) {
            return Err(PreparedError::InvalidFasta(
                "duplicate FASTA identifier".into(),
            ));
        }
        dictionary_output.write_all(name).map_err(index_io)?;
        writeln!(dictionary_output, "\t{}", record.length()).map_err(index_io)?;
        total_bases = total_bases
            .checked_add(record.length())
            .ok_or_else(|| limit("base count overflow"))?;
        sequence_count += 1;
        // Upstream provides only write_index; a single-record index keeps memory
        // independent of total FASTA size and avoids a BioV FAI serializer.
        fai::io::Writer::new(&mut fai_output)
            .write_index(&fai::Index::from(vec![record]))
            .map_err(index_io)?;
    }
    if sequence_count == 0 {
        return Err(PreparedError::InvalidFasta(
            "FASTA has no sequence records".into(),
        ));
    }
    fai_output.flush().map_err(index_io)?;
    dictionary_output.flush().map_err(index_io)?;
    Ok(Indexed {
        sequence_count,
        total_bases,
        source_sha256: format!("{:x}", source.digest.finalize()),
        source_bytes: source.bytes,
        fai: fai_output.identity("sequences.fai"),
        dictionary: dictionary_output.identity("sequences.tsv"),
    })
}
fn index_io(error: io::Error) -> PreparedError {
    match error.kind() {
        io::ErrorKind::FileTooLarge => limit(&error.to_string()),
        io::ErrorKind::InvalidInput | io::ErrorKind::InvalidData | io::ErrorKind::UnexpectedEof => {
            PreparedError::InvalidFasta(error.to_string())
        }
        _ => crate::io("stream FASTA or prepared output", error),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn upstream_output_is_independent_of_input_buffer_boundaries() {
        for source in [
            b">sq0 description\nACGT\nAC\n>sq1\nNNNN\n".as_slice(),
            b">sq0\r\nACGT\r\nAC\r\n>sq1\r\nNNNN",
            b">sq0\nACGT",
        ] {
            let mut expected = None;
            for capacity in [1, 2, 3, 4, 8, 64, IO_BUFFER_BYTES] {
                let mut fai = Vec::new();
                let mut dict = Vec::new();
                let actual = index_with_capacity(source, &mut fai, &mut dict, capacity).unwrap();
                if let Some(expected) = &expected {
                    assert_eq!(&(actual, fai, dict), expected);
                } else {
                    expected = Some((actual, fai, dict));
                }
            }
        }
    }
    #[test]
    fn malformed_or_unsupported_inputs_reject_at_every_buffer_capacity() {
        for source in [
            b">sq0\nAC\rGT\n".as_slice(),
            b">sq0\nAC>GT\n",
            b">sq0\nAC GT\n",
            b">sq0\nAC1GT\n",
            b">sq0\nAC-GT\n",
            b">sq0\nACUGT\n",
            b">sq0\nAC\0GT\n",
            b">sq0\nACGT\n\n",
            b">sq0\n",
            b"",
            b"\n",
            b">sq0\nAC\n>sq0\nGT\n",
            b">sq0\nAC\nACGT\n",
            b">sq0\nACGT\nA\nACGT\n",
            b"> sq0\nACGT\n",
            b">sq0\nACGT\r",
            b">sq0\nACGT\nACGT\r\nACGT\n",
            b">\xff\nACGT\n",
            b"\x1f\x8bcompressed",
        ] {
            for capacity in [1, 2, 3, 8, 64] {
                assert!(
                    index_with_capacity(source, io::sink(), io::sink(), capacity).is_err(),
                    "accepted {source:?} capacity {capacity}"
                );
            }
        }
    }
    #[test]
    fn headers_sequence_lines_and_identifiers_are_bounded() {
        let mut source = vec![b'>'; MAX_LINE_BYTES + 1];
        source.extend_from_slice(b"\nACGT\n");
        assert!(matches!(
            index(&source[..], io::sink(), io::sink()),
            Err(PreparedError::Limit(_))
        ));
        let mut source = b">sq0\n".to_vec();
        source.extend(std::iter::repeat_n(b'A', MAX_LINE_BYTES + 1));
        assert!(matches!(
            index(&source[..], io::sink(), io::sink()),
            Err(PreparedError::Limit(_))
        ));
        let mut source = b">".to_vec();
        source.extend(std::iter::repeat_n(b'A', MAX_ID_BYTES + 1));
        source.extend_from_slice(b"\nACGT\n");
        assert!(matches!(
            index(&source[..], io::sink(), io::sink()),
            Err(PreparedError::Limit(_))
        ));
    }
    #[test]
    fn identifiers_and_case_are_preserved_with_standard_fai_serialization() {
        let source = b">0001:alt description\tmore\nacgTRYSW\nKMBDHVN\n>second\nNN\n";
        let mut fai = Vec::new();
        let mut dict = Vec::new();
        let result = index(&source[..], &mut fai, &mut dict).unwrap();
        assert_eq!(result.sequence_count, 2);
        assert_eq!(result.total_bases, 17);
        assert_eq!(dict, b"sequence_id\tlength\n0001:alt\t15\nsecond\t2\n");
        let parsed = fai::io::Reader::new(&fai[..]).read_index().unwrap();
        assert_eq!(parsed.as_ref()[0].name(), b"0001:alt");
        assert_eq!(parsed.as_ref()[0].length(), 15);
        assert_eq!(parsed.as_ref()[0].position(), 27);
    }
    #[test]
    fn record_count_identifier_charge_and_output_bytes_have_independent_limits() {
        let mut source = Vec::new();
        for i in 0..=MAX_RECORDS {
            writeln!(source, ">{i}\nA").unwrap();
        }
        assert!(
            matches!(index(&source[..], io::sink(), io::sink()), Err(PreparedError::Limit(message)) if message.contains("100000"))
        );
        let mut source = Vec::new();
        for i in 0..40_000 {
            writeln!(source, ">{i:06}{}\nA", "x".repeat(300)).unwrap();
        }
        assert!(
            matches!(index(&source[..], io::sink(), io::sink()), Err(PreparedError::Limit(message)) if message.contains("metadata charge"))
        );
        let mut output = HashWriter::new(io::sink());
        let chunk = [0u8; IO_BUFFER_BYTES];
        for _ in 0..MAX_METADATA_BYTES / IO_BUFFER_BYTES {
            output.write_all(&chunk).unwrap();
        }
        assert_eq!(output.bytes, MAX_METADATA_BYTES as u64);
        assert_eq!(
            output.write(b"x").unwrap_err().kind(),
            io::ErrorKind::FileTooLarge
        );
    }
    #[test]
    fn exact_physical_line_bound_accepts_a_long_unwrapped_final_sequence() {
        let mut source = b">a\n".to_vec();
        source.extend(std::iter::repeat_n(b'A', MAX_LINE_BYTES));
        let result = index(&source[..], io::sink(), io::sink()).unwrap();
        assert_eq!(result.total_bases, MAX_LINE_BYTES as u64);
        source.push(b'\n');
        assert!(matches!(
            index(&source[..], io::sink(), io::sink()),
            Err(PreparedError::Limit(_))
        ));
    }
}
