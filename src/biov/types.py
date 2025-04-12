from typing import Annotated

from Bio.Seq import Seq as _Seq
from pydantic import PlainSerializer, PlainValidator, WithJsonSchema

Seq = Annotated[
    _Seq,
    PlainValidator(_Seq),
    PlainSerializer(str),
    WithJsonSchema({"type": "string"}, mode="serialization"),
]
