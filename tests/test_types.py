from pydantic import BaseModel, TypeAdapter

from biov.types import Seq


class A(BaseModel):
    seq: Seq


def test_seq_validation_and_serialization():
    acgt = TypeAdapter(Seq).validate_python("ACGT")
    assert acgt == "ACGT"
    assert acgt.reverse_complement() == "ACGT"
    a = A.model_validate({"seq": "ACGT"})
    assert a.seq == acgt
    assert a.model_dump() == {"seq": "ACGT"}
    assert a.model_dump_json() == '{"seq":"ACGT"}'
