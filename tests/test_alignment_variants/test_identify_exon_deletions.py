import pytest
from kir_annotator.alignment_variants import identify_exon_deletions 
from collections import namedtuple

# Mock objects and data for testing
class MockRegion:
    def __init__(self, name, start, end, is_exon, pseudo, seq):
        self.name = name
        self.start = start
        self.end = end
        self.is_exon = is_exon
        self.pseudo = pseudo
        self.seq = seq

gene = namedtuple('Gene', ['name'])

class MockAllele:
    def __init__(self, regions):
        self.regions = regions
        self.gene = gene('KIR3DL1')
    
    def region(self, i):
        """Return the region name and offset of gene location."""

        st = 0
        for r in self.regions.values():
            if st <= i < st + len(r.seq):
                return r, i - st
            st += len(r.seq)
        # assert False, i
        return None, None


# Actual tests
@pytest.mark.parametrize("pos, op, expected_results", [
    # Test case 1: Deletion within a single exon
    (0, 'delAAA', [(0, 'delAAA')]),

    # Test case 2: Deletion spanning multiple exons
    (18, 'delTTGGGGCCCCCCCCCC', [(18, 'delTT'), (32, 'delCC')]),

    # Test case 3: Deletion outside exons (should return empty)
    (20, 'delGG', []),

    # Test case 4: Deletion of all regions
    (0, 'delAAAAAAAAAATTTTTTTTTTGGGGCCCCCCCCCCCCGGGGTTTTAAAA', [(0, 'delAAAAAAAAAATTTTTTTTTT'), (32, 'delCCCCGGGGTTTTAAAA')]),
])
def test_identify_exon_deletions(pos, op, expected_results):
    # Setup mock data for wildtype allele
    mock_regions = [
        MockRegion('e1', 0, 20, True, False, 'AAAAAAAAAATTTTTTTTTT'),  # Exon
        MockRegion('i1', 20, 32, False, False, 'GGGGCCCCCCCC'),       # Intron
        MockRegion('e2', 32, 48, True, False, 'CCCCGGGGTTTTAAAA'),    # Exon
        # Add more mock regions as needed for your tests
    ]
    mock_allele = MockAllele({1: mock_regions[0], 2: mock_regions[1], 3: mock_regions[2]})  # Adjust based on your structure

    # Run the function with the test case
    result = identify_exon_deletions(pos, op, mock_allele)

    # Check that the result matches the expected output
    assert result == expected_results, f"Failed for pos={pos}, op={op}"
