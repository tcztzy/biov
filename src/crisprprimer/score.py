"""Deterministic existing scores with the original bundled lookup tables."""

from pathlib import Path
from typing import cast

import pandas as pd

BASE_PATH = Path(__file__).parent
mm_scores = pd.read_table(
    BASE_PATH / "cfd.mm.scores.cas9.txt", sep=" ", names=["MM", "score"], index_col=0
)["score"]
pam_scores = pd.read_table(
    BASE_PATH / "cfd.pam.scores.cas9.txt", sep=" ", names=["PAM", "score"], index_col=0
)["score"]


def cfd_score(spacer: str, protospacer: str, pam: str) -> float:
    """Calculate the existing CFD score from the bundled lookup tables.

    Returns:
        The existing computed CFD score, not measured specificity.

    Raises:
        ValueError: If spacer or PAM lengths do not match the score model.
    """
    score = 1.0
    if len(spacer) != 20 or len(protospacer) != 20 or len(pam) != 2:
        raise ValueError("CFD requires 20-base spacers and a 2-base PAM suffix")
    for i, (s, p) in enumerate(zip(spacer, protospacer)):
        if s != p:
            score *= cast(float, mm_scores[f"{s}{p}{i + 1}"])
    score *= cast(float, pam_scores.get(pam, 0))
    return score


def distance_to_cds_start_score(cds: pd.DataFrame, start: int, end: int) -> float:
    """Calculate the score of a CDS region.

    Returns:
        The existing distance score between 0.5 and 1.
    """
    # Cover all CDS
    if start <= cds["start"].min() and end >= cds["end"].max():
        return 1.0
    # part of design is outside
    if (
        start < cds["end"].max()
        and end > cds["start"].min()
        and not (start > cds["start"].min() and end < cds["end"].max())
    ):
        return 0.5
    olen = cds["end"].max() - cds["start"].min()
    cent = (start + end) / 2
    dist = (
        abs(cent - cds["start"].min())
        if (cds["strand"] == "+").all()
        else abs(cent - cds["end"].max())
    )
    if olen <= 10:
        return 1.0
    if dist <= 10:
        return 1.0
    return max(0.5, min(1 - dist / olen, 1.0))


def length_score(start: int, end: int, dist: int = 3) -> float:
    """Calculate a score based on sequence length relative to target length.

    Args:
        start: Start position
        end: End position
        dist: Distance parameter (default=3)

    Returns:
        float: Score between 0.0 and 1.0
    """
    length = end - start
    target_length = 60 + (dist * 2)  # 60 is the window size between two cutting sites
    if length < 4:
        return 0.0
    deviation = abs(length - target_length)
    return max(0.5, 1 - (deviation / 30))


def pair_score(row: pd.Series, cds: pd.DataFrame, dist: int) -> float:
    """Combine the existing CDS-distance and interval-length scores.

    Returns:
        The product of the two existing component scores.
    """
    start = cast(int, row["start"])
    end = cast(int, row["end"])
    return distance_to_cds_start_score(cds, start, end) * length_score(start, end, dist)
