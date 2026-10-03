"""Typed native batch sequence boundary."""

class SequenceValidationError(ValueError):
    """A string contains symbols outside its declared alphabet."""

def normalize_sequences(values: list[str | None], *, kind: str) -> list[str | None]:
    """Normalize a nullable sequence batch with an explicit kind."""

def reverse_complements(values: list[str | None], *, kind: str) -> list[str | None]:
    """Reverse-complement a nullable DNA or RNA batch."""

def sequence_lengths(values: list[str | None], *, kind: str) -> list[int | None]:
    """Count validated symbols in nullable DNA, RNA or protein sequences."""

def weighted_gc_fractions(values: list[str | None], *, kind: str) -> list[float | None]:
    """Compute the mean GC probability of IUPAC nucleotide codes."""
