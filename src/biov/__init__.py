"""Next-generation development experience for computational molecular biology."""

from . import _patch  # noqa

from .config import settings as settings
from .dataframe import BioDataFrame as BioDataFrame
from .io.fastx import read_fasta as read_fasta
from .io.gff import read_gff3 as read_gff3
from .seq import Seq as Seq
from .seq import SequenceArray as SequenceArray
from .seq import SequenceDtype as SequenceDtype
from .seq import SequenceValidationError as SequenceValidationError

__all__ = [
    "BioDataFrame",
    "Seq",
    "SequenceArray",
    "SequenceDtype",
    "SequenceValidationError",
    "read_fasta",
    "read_gff3",
    "settings",
]
