from __future__ import annotations

import itertools
from typing import Callable
from urllib.parse import unquote

import pandas as pd
from pandas import DataFrame
from pandas._typing import FilePath, ReadCsvBuffer

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

    def attributes_to_columns(self) -> GFFDataFrame:
        """Saving each attribute-tag to a single column.

        Attribute column will be split by the tags in the single columns.
        For this method only a pandas DataFrame and not a Gff3DataFrame
        will be returned. Therefore, this data frame can not be saved as
        gff3 file.
        """
        df = self.copy()[GFF3_COLUMNS]
        attributes_dict = self["attributes"].apply(
            lambda attributes: dict(kv.split("=") for kv in attributes.split(";"))
        )
        all_attributes = set(
            itertools.chain.from_iterable(attributes_dict.apply(lambda d: d.keys()))
        )
        attributes = [
            a
            for a in itertools.chain(
                GFF3_ATTRIBUTES, sorted(all_attributes - set(GFF3_ATTRIBUTES))
            )
            if a in all_attributes
        ]
        for attribute in attributes:
            df[attribute] = attributes_dict.apply(lambda d: d.get(attribute)).apply(
                lambda v: unquote(v) if isinstance(v, str) else v
            )
        return df  # pyright: ignore


def read_gff3(input_file: FilePath | ReadCsvBuffer[bytes] | ReadCsvBuffer[str]):
    return GFFDataFrame(pd.read_table(input_file, comment="#", names=GFF3_COLUMNS))
