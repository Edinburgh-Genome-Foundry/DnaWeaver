import dnaweaver.biotools as bt
import pytest


def test_largest_common_substring():
    seqA = "-----oooooooo"
    seqB = "oooooo-----tttt"

    assert bt.largest_common_substring(seqA, seqA, 80) == (0, 12)
    assert bt.largest_common_substring(seqA, seqB, 80) == (5, 11)
    assert bt.largest_common_substring(seqA, seqB, 5) == (5, 11)
    assert bt.largest_common_substring(seqA, seqB, 4) is False


@pytest.mark.parametrize("seq, expected", [("AATCG", 0.4), ("aatgc", 0.4)])
def test_gc_content(seq, expected):
    assert bt.gc_content(seq) == expected
