"""Read FASTA and FASTQ file."""

import fsspec
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pandas._typing import FilePath, ReadBuffer

from ._preprocess import preprocessing


def read_fasta(
    filepath_or_buffer: FilePath | ReadBuffer,
) -> dict[str, SeqRecord]:
    """Read FASTA sequence.

    Returns:
        a dict key as sequence id and value as SeqRecord.
    """
    filepath_or_buffer, _ = preprocessing(filepath_or_buffer)
    if isinstance(filepath_or_buffer, str):
        with fsspec.open(filepath_or_buffer, "rt", compression="infer") as b:
            return SeqIO.to_dict(SeqIO.parse(b, "fasta"))
    else:
        return SeqIO.to_dict(SeqIO.parse(filepath_or_buffer, "fasta"))
