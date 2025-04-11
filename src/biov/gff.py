from __future__ import annotations

from typing import Callable
from urllib.parse import unquote

import pandas as pd
from pandas import DataFrame
from pandas._typing import FilePath, ReadCsvBuffer

from .config import settings

GFF3_COLUMNS = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]
GFF3_ATTRIBUTES = (
    "ID",
    "Name",
    "Alias",
    "Parent",
    "Target",
    "Gap",
    "Derives_from",
    "Note",
    "Dbxref",
    "Ontology_term",
    "Is_circular",
)


class GFFDataFrame(DataFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for col in GFF3_COLUMNS:
            if col not in self.columns:
                raise AttributeError(f"Column '{col}' is required")

    @property
    def _constructor(self) -> Callable[..., GFFDataFrame]:
        return GFFDataFrame

    def to_gff3(self, gff_file: FilePath) -> None:
        df = self[GFF3_COLUMNS]
        gff_feature = df.to_csv(sep="\t", index=False, header=False)
        with open(gff_file, "w") as fh:
            fh.write("##gff-version 3\n")
            fh.write(gff_feature)

    @property
    def attributes(self) -> DataFrame:
        df = pd.json_normalize(
            self["attributes"].apply(
                lambda attributes: dict(kv.split("=") for kv in attributes.split(";"))
            )  # type: ignore
        ).map(lambda v: unquote(v) if isinstance(v, str) else v)
        df.index = self.index
        order: dict[str, int] = dict((a, i) for i, a in enumerate(GFF3_ATTRIBUTES))
        columns = sorted(list(df.columns), key=lambda c: order.get(c, len(c)))
        return df[columns]  # pyright: ignore


def read_gff3(
    input_file: FilePath | ReadCsvBuffer[bytes] | ReadCsvBuffer[str], **kwargs
):
    """Read GFF3 files.

    Parameters
    ----------
    input_file : str | os.PathLike | ReadableBuffer
        support fsspec chain (available in pandas 3.0)
    **kwargs
        will pass to `pd.read_table`

    Returns
    -------
    GFFDataFrame

    Raises
    ------
    ValueError
        if provided 'comment' parameter
    """
    for param in ("comment", "na_values"):
        if param in kwargs:
            raise ValueError(f"Parameter '{param}' is not allowed")
    kwargs["comment"] = "#"
    kwargs["na_values"] = "."
    if "names" not in kwargs:
        kwargs["names"] = GFF3_COLUMNS
    if isinstance(input_file, str):
        *protocols, path = input_file.split("::")
        # https://github.com/pandas-dev/pandas/pull/60100
        if any([protocol.startswith("tar://") for protocol in protocols]):
            kwargs["compression"] = None
        if (
            settings.cache_http
            and path.startswith(("https://", "http://"))
            and len(protocols) == 0
            or "filecache" != protocols[-1]
        ):
            input_file = "::".join([*protocols, "filecache", path])
    return GFFDataFrame(pd.read_table(input_file, **kwargs))
