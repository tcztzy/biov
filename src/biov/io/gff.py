"""Generic Feature Format Version 3

Specification:
    https://github.com/the-sequence-ontology/specifications/blob/master/gff3.md
"""

import os
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal
from urllib.parse import unquote

import pandas as pd
from pandas._typing import FilePath, ReadCsvBuffer, WriteBuffer

from ._preprocess import preprocessing

if TYPE_CHECKING:
    from ..dataframe import BioDataFrame


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


GFF_COLUMNS = [
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


def read_gff3(
    filepath_or_buffer: FilePath | ReadCsvBuffer[bytes] | ReadCsvBuffer[str],
    explode_attributes: bool = True,
    **kwargs: Any,
) -> "BioDataFrame":
    """Read GFF3 files.

    Args:
        filepath_or_buffer: support fsspec chain (available in pandas 3.0)
        explode_attributes: explode attributes into columns, notice that attributes with same name as primary columns (such as `score`) will be ignored.
        **kwargs: will pass to `pd.read_table`

    Returns:
        A BioDataFrame with at least 9 columns

    Examples:
        You can directly pass HTTP url.

        >>> read_gff3("https://rice.uga.edu/osa1r7_download/osa1_r7.asm.repeat_masked.gff3.gz")

        You can also pass fsspec chain

        >>> read_gff3("tar://IRGSP-1.0_representative/transcripts.gff::https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_2025-03-19.tar.gz")
    """

    for param in ("comment", "na_values"):
        if param in kwargs:
            raise ValueError(f"Parameter '{param}' is not allowed")
    kwargs["comment"] = "#"
    kwargs["na_values"] = "."
    if "names" not in kwargs:
        kwargs["names"] = GFF_COLUMNS
    names: list[str] = list(kwargs["names"])
    filepath_or_buffer, kwargs = preprocessing(filepath_or_buffer, **kwargs)

    def explode(attributes: pd.Series):
        df = pd.json_normalize(
            attributes.apply(
                lambda attrs: dict(
                    [unquote(i) for i in kv.split("=", maxsplit=1)]
                    for kv in attrs.split(";")
                )
                if isinstance(attrs, str)
                else {}
            )  # type: ignore
        )
        return df[[c for c in df.columns if c not in names]]

    df: pd.DataFrame = pd.read_table(filepath_or_buffer, dtype={}, **kwargs)
    df["start"] -= 1
    if explode_attributes:
        df = pd.concat(
            [df[names[:-1]], explode(df[names[-1]])],  # pyright: ignore
            axis=1,
        )
    from ..dataframe import BioDataFrame

    bdf = BioDataFrame(df)
    bdf._gff_columns = names
    return bdf


def to_gff3(
    self: "BioDataFrame",
    path: FilePath | WriteBuffer[bytes] | None = None,
    *,
    attributes: str | list[str] | Literal[True] = True,
) -> str | None:
    names = self._gff_columns
    if len(names) != 9:
        raise ValueError("Mal-formatted BioDataFrame for export to GFF3.")
    if attributes is True:
        attributes = [c for c in self.columns if c not in names]
    elif isinstance(attributes, str):
        attributes = [attributes]
    df = self.copy()
    df[names[3]] += 1
    defaults = {
        1: "biov",  # source
        2: "gene",  # type
        5: pd.NA,  # score
        6: pd.NA,  # strand
        7: pd.NA,  # phase
    }
    for c in range(8):  # except attributes
        if (name := names[c]) not in self.columns:
            if (default := defaults.get(c)) is not None:
                df[name] = default
            else:
                raise ValueError(f"Column `{name}` is required for exporting to GFF.")
    df = df[names[:-1]]
    if len(attributes) == 0:
        df["attributes"] = pd.NA
    else:
        df["attributes"] = [
            ";".join(
                f"{quote(k)}={quote(str(v))}"
                for k, v in r.items()
                if isinstance(k, str) and not pd.isna(v)
            )
            for r in self[attributes].to_dict(orient="records")  # pyright: ignore
        ]
    gff_feature = df.to_csv(sep="\t", na_rep=".", index=False, header=False)
    content = f"""##gff-version 3\n{gff_feature}"""
    if path is None:
        return content
    if isinstance(path, str | os.PathLike):
        Path(path).write_text(content)
        return None
    path.write(content.encode())
    return None


class GFFMixin:
    _gff_columns: list[str] = GFF_COLUMNS

    to_gff3 = to_gff3
    to_gff = to_gff3
