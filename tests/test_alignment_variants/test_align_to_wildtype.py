import pytest
from kir_annotator.alignment_variants import align_to_wildtype

# Test cases
@pytest.mark.parametrize("assembly_seq, wildtype_seq, expected", [
    ("AAGT", "ACGT", ('AAGT', [(1, '='), (1, 'X'), (2, '=')], False)),
    ("AACGT", "ACGT", ('AACGT', [(1, 'D'), (4, '=')], False)),
    ("ACGTT", "ACGT", ('ACGTT', [(3, '='), (1, 'D'), (1, '=')], False)),
    ("ACGT", "AGT", ('ACGT', [(1, '='), (1, 'D'), (2, '=')], False)),
    ("AGT", "ACGT", ('AGT', [(1, 'I'), (1, 'X'), (2, '=')], False)),
    ("ACAGT", "ACGT", ('ACAGT', [(2, '='), (1, 'D'), (2, '=')], False))
])
def test_align_to_wildtype(assembly_seq, wildtype_seq, expected):
    """
    Test align_to_wildtype function with various sequences.

    Args:
        assembly_seq (str): Assembly sequence.
        wildtype_seq (str): Wildtype sequence.
        expected (Tuple): Expected result.
    """
    assert align_to_wildtype(assembly_seq, wildtype_seq) == expected


def test_empty_sequences():
    """
    Test align_to_wildtype function with empty sequences.
    """
    with pytest.raises(ValueError):
        align_to_wildtype("", "ACGT")
    with pytest.raises(ValueError):
        align_to_wildtype("ACGT", "")