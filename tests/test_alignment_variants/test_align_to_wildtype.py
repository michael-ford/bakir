import pytest
from kir_annotator.alignment_variants import align_to_wildtype

# Test cases
@pytest.mark.parametrize("assembly_seq, wildtype_seq, is_reversed, expected", [
    ("AAGT", "ACGT", False, ('AAGT', [(1, '='), (1, 'X'), (2, '=')])),
    ("AACGT", "ACGT", False, ('AACGT', [(1, 'D'), (4, '=')])),
    ("ACGTT", "ACGT", False, ('ACGTT', [(3, '='), (1, 'D'), (1, '=')])),
    ("ACGT", "AGT", False, ('ACGT', [(1, '='), (1, 'D'), (2, '=')])),
    ("AGT", "ACGT", False, ('AGT', [(1, 'I'), (1, 'X'), (2, '=')])),
    ("ACAGT", "ACGT", False, ('ACAGT', [(2, '='), (1, 'D'), (2, '=')]))
])
def test_align_to_wildtype(assembly_seq, wildtype_seq, is_reversed, expected):
    """
    Test align_to_wildtype function with various sequences.

    Args:
        assembly_seq (str): Assembly sequence.
        wildtype_seq (str): Wildtype sequence.
        expected (Tuple): Expected result.
    """
    assert align_to_wildtype(assembly_seq, wildtype_seq, is_reversed) == expected


def test_empty_sequences():
    """
    Test align_to_wildtype function with empty sequences.
    """
    with pytest.raises(ValueError):
        align_to_wildtype("", "ACGT", False)
    with pytest.raises(ValueError):
        align_to_wildtype("ACGT", "", False)