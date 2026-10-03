"""Rice identifier lookup using the published mapping and bundled data."""

import re
from typing import Literal

import fsspec

RAP2MSU: dict[str, list[str] | None] = {}
MSU2RAP: dict[str, str | None] = {}
MAP_FILE = "https://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU_2025-03-19.txt.gz"
_mapping_loaded = False


def _load_mapping() -> None:
    global _mapping_loaded
    if _mapping_loaded:
        return

    with fsspec.open(
        f"filecache::{MAP_FILE}",
        mode="rt",
        compression="gzip",
    ) as f:
        for line in f:
            rap_id, msu_ids = line.strip().split("\t")
            if rap_id == "None":
                for i in msu_ids.split(","):
                    MSU2RAP[i] = None
            elif msu_ids == "None":
                RAP2MSU[rap_id] = None
            else:
                RAP2MSU[rap_id] = ids = msu_ids.split(",")
                for i in ids:
                    MSU2RAP[i] = rap_id

    _mapping_loaded = True


def guess_id_system(id: str) -> Literal["MSU", "RAP"] | None:
    """Recognize a rice RAP or MSU identifier.

    Returns:
        RAP, MSU, or None for an unrecognized identifier.
    """
    if re.match(r"Os(?:0[1-9]|1[012])g\d{7}", id):
        return "RAP"
    elif re.match(r"LOC_Os(?:0[1-9]|1[012])g\d{5}(\.\d+)?", id):
        return "MSU"
    else:
        return None


def convert(id: str) -> str | list[str] | None:
    """Look up the published RAP-MSU mapping.

    Returns:
        Mapped identifiers, or None if no mapping exists.

    Raises:
        ValueError: If the input is neither a RAP nor MSU identifier.
    """
    id_system = guess_id_system(id)
    if id_system is None:
        raise ValueError(f"Your input id {id} is not RAP id or MSU id.")
    _load_mapping()
    if id_system == "RAP":
        return RAP2MSU.get(id)
    else:
        if "." not in id:
            id += ".1"
        return MSU2RAP.get(id)
