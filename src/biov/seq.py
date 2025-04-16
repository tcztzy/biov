"""Types for BioV."""

from typing import Any

from Bio.Seq import Seq as _Seq
from pydantic import (
    GetCoreSchemaHandler,
)
from pydantic_core import CoreSchema, core_schema


class Seq(_Seq):
    """Extended from Bio.Seq.Seq for pydantic validation and serialization."""

    @classmethod
    def __get_pydantic_core_schema__(
        cls, source_type: Any, handler: GetCoreSchemaHandler
    ) -> CoreSchema:
        """Validation and serialization for Seq.

        Returns:
            validation and serialization schema
        """
        return core_schema.no_info_plain_validator_function(
            _Seq,
            serialization=core_schema.to_string_ser_schema(),
        )
