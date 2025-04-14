from typing import Any

from Bio.Seq import Seq as _Seq
from pydantic import (
    GetCoreSchemaHandler,
)
from pydantic_core import CoreSchema, core_schema


class Seq(_Seq):
    @classmethod
    def __get_pydantic_core_schema__(
        cls, source_type: Any, handler: GetCoreSchemaHandler
    ) -> CoreSchema:
        return core_schema.no_info_plain_validator_function(
            _Seq,
            serialization=core_schema.plain_serializer_function_ser_schema(
                str, return_schema=core_schema.str_schema()
            ),
        )
