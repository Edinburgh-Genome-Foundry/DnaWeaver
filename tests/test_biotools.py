import dnaweaver.biotools as bt


def test_largest_common_substring():

    seqA = "-----oooooooo"
    seqB = "oooooo-----tttt"

    assert bt.largest_common_substring(seqA, seqA, 80) == (0, 12)
    assert bt.largest_common_substring(seqA, seqB, 80) == (5, 11)
    assert bt.largest_common_substring(seqA, seqB, 5) == (5, 11)
    assert bt.largest_common_substring(seqA, seqB, 4) is False
