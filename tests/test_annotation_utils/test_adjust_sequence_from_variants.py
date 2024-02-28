import pytest
from kir_annotator.annotation_utils import adjust_sequence_from_variants

# Mock classes
class MockWildType:
    def __init__(self, seq):
        self.seq = seq

class MockGene:
    def __init__(self, wildtype, mutations):
        self.wildtype = wildtype
        self.mutations = mutations

wildtype_example = MockWildType("ATCGATCG")
common_mutations = [(4, "insGG"), (6, "delCC")]  # These mutations are in the middle and should not affect start/end adjustments

# Test cases
@pytest.mark.parametrize("variants, expected", [
    ([(0, "insGG")], (2, 0)),  # 2 bp insertion at the beginning
    ([(0, "delAT")], (-2, 0)), # 2 bp deletion at the beginning
    ([(8, "insTT")], (0, -2)),  # 2 bp insertion at the end
    ([(6, "delCG")], (0, 2)), # 2 bp deletion at the end
    ([(0, "insGG"), (0, "insAA")], pytest.raises(ValueError)),  # Multiple start indels
    ([(8, "delCG"), (6, "delCG")], pytest.raises(ValueError)),  # Multiple end indels
    ([], (0, 0)),  # No indels
    ([(3, "insGG")], (0, 0)),  # Indel in the middle
])
def test_adjust_sequence_from_variants(variants, expected):
    gene_example = MockGene(wildtype_example, common_mutations)
    if isinstance(expected, tuple):
        assert adjust_sequence_from_variants(variants, gene_example) == expected
    else:
        with expected:
            adjust_sequence_from_variants(variants, gene_example)
