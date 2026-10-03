"""Existing nuclease matching through Biopython and validated configuration."""

from enum import Enum
from typing import Any

import pandas as pd
from Bio.Seq import MutableSeq
from Bio.Seq import Seq as BioSeq
from Bio.SeqUtils import nt_search
from pydantic import BaseModel, field_validator

from biov import BioDataFrame, Seq


class Nuclease(BaseModel, frozen=True):
    """Nuclease model."""

    pam: Seq
    cut_sites: tuple[int, int]
    """PAM sequence"""
    spacer_range: tuple[int, int]

    @field_validator("cut_sites", mode="before")
    @classmethod
    def validate_cut_sites(cls, value: Any) -> Any:
        """Expand one cut position to the two-strand representation.

        Returns:
            The value to validate as a pair.
        """
        if isinstance(value, int):
            return (value, value)
        return value

    @field_validator("pam")
    @classmethod
    def validate_pam(cls, value: Seq) -> Seq:
        """Reject an empty PAM.

        Returns:
            The nonempty PAM.

        Raises:
            ValueError: If the PAM is empty.
        """
        if len(value) <= 0:
            raise ValueError("Zero length PAM is not allowed.")
        return value

    @field_validator("spacer_range")
    @classmethod
    def validate_spacer_range(cls, value: tuple[int, int]) -> tuple[int, int]:
        """Require a positive spacer interval length.

        Returns:
            The increasing interval.

        Raises:
            ValueError: If the interval is not increasing.
        """
        if value[1] <= value[0]:
            raise ValueError("spacer_range should be (start, end) and start < end")
        return value

    @property
    def pam_side(self) -> str:
        """Identify which side of the spacer contains the PAM."""
        return "3prime" if self.spacer_range[1] <= 0 else "5prime"

    @property
    def spacer_length(self) -> int:
        """The configured spacer interval length."""
        return self.spacer_range[1] - self.spacer_range[0]

    @property
    def prototype(self) -> str:
        """Build the existing ambiguous sequence-search template."""
        proto = str(self.pam)
        if self.spacer_range[0] < 0:
            proto = "N" * abs(self.spacer_range[0]) + proto
        if self.spacer_range[1] > len(self.pam):
            proto = proto + "N" * (self.spacer_range[1] - len(self.pam))
        return proto

    def find_spacers_on(self, seq: BioSeq | MutableSeq) -> BioDataFrame:
        """Find both-strand matches with Biopython ambiguity matching.

        Returns:
            Matches in 0-based, end-exclusive coordinates.
        """
        start = pd.Series(nt_search(str(seq), self.prototype)[1:])
        df = BioDataFrame()
        if len(start) > 0:
            df["start"] = start + (len(self.pam) if self.pam_side == "5prime" else 0)
            df["end"] = df["start"] + self.spacer_length
            df["pam"] = (start if self.pam_side == "5prime" else df["end"]).apply(
                lambda s: str(seq[s : s + len(self.pam)])
            )
            df["strand"] = "+"
            df["spacer"] = df.apply(
                lambda row: str(seq[row["start"] : row["end"]]),
                axis=1,
            )

        start = pd.Series(
            nt_search(str(seq), str(Seq(self.prototype).reverse_complement()))[1:]
        )
        reverse_df = BioDataFrame()
        if len(start) > 0:
            reverse_df["start"] = start + (
                len(self.pam) if self.pam_side == "3prime" else 0
            )
            reverse_df["end"] = reverse_df["start"] + self.spacer_length
            reverse_df["pam"] = (
                start if self.pam_side == "3prime" else reverse_df["end"]
            ).apply(lambda s: str(seq[s : s + len(self.pam)].reverse_complement()))
            reverse_df["strand"] = "-"
            reverse_df["spacer"] = reverse_df.apply(
                lambda row: str(seq[row["start"] : row["end"]].reverse_complement()),
                axis=1,
            )
        if df.empty and reverse_df.empty:
            return df
        df = BioDataFrame(
            pd.concat([d for d in [df, reverse_df] if not d.empty], ignore_index=True)
        )
        return df


SpCas9 = Nuclease(cut_sites=(-3, -3), pam=Seq("NGG"), spacer_range=(-20, 0))
AsCas12 = Nuclease(cut_sites=(22, 27), pam=Seq("TTTV"), spacer_range=(4, 27))


class Nucleases(Enum):
    """Existing named nuclease configurations."""

    SpCas9 = "SpCas9"
    AsCas12a = "AsCas12a"

    @property
    def model(self) -> Nuclease:
        """The selected nuclease configuration."""
        return SpCas9 if self is Nucleases.SpCas9 else AsCas12
