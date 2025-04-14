from __future__ import annotations

from typing import Callable, Literal
from urllib.parse import unquote

import pandas as pd
from pandas import DataFrame
from pandas._typing import FilePath, ReadCsvBuffer

from .config import settings


def quote(s: str) -> str:
    return s.translate(
        str.maketrans(
            {  # type: ignore[arg-type]
                ";": "%3B",
                "=": "%3D",
                "&": "%26",
                ",": "%2C",
                chr(0x7F): "%7F",
                **{chr(i): f"%{i:x}".upper() for i in range(0x1F)},
            }
        )
    )


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

    def to_gff3(
        self,
        gff_file: FilePath | None = None,
        *,
        extra: Literal["ignore", "merge", "inplace"] = "ignore",
    ) -> str | None:
        df = self[GFF3_COLUMNS].copy()
        if extra != "ignore":
            extra_fields = [
                {k: v for k, v in e.items() if k not in GFF3_COLUMNS}
                for e in self.to_dict(orient="records")
            ]
            if extra == "merge":
                attributes_field = [
                    a | e
                    for a, e in zip(
                        self.attributes.to_dict(orient="records"), extra_fields
                    )
                ]
            else:
                attributes_field = extra_fields
            df["attributes"] = [
                ";".join(
                    f"{quote(k)}={quote(str(v))}"
                    for k, v in r.items()
                    if isinstance(k, str) and not pd.isna(v)
                )
                for r in attributes_field
            ]
        gff_feature = df.to_csv(sep="\t", index=False, header=False)
        if gff_file is None:
            return f"""##gff-version 3\n{gff_feature}"""
        with open(gff_file, "w") as fh:
            fh.write("##gff-version 3\n")
            fh.write(gff_feature)
            return None

    @property
    def version(self) -> str:
        return "3"

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

    def attributes_to_columns(self) -> GFFDataFrame:
        attributes = self.attributes
        df = self.copy()
        for c in attributes.columns:
            if c not in GFF3_COLUMNS:
                df[c] = attributes[c]
        return df


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
