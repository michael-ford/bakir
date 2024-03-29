import pytest
from collections import namedtuple
from kir_annotator.alignment_variants import identify_functional_deletion

# Define a simple structure for regions as before
Region = namedtuple('Region', ['name', 'is_exon', 'pseudo'])

# Fixtures for different types of regions
@pytest.fixture
def functional_exon():
    return Region(name="exon1", is_exon=True, pseudo=False)

@pytest.fixture
def pseudo_exon():
    return Region(name="pseudoexon1", is_exon=True, pseudo=True)

@pytest.fixture
def non_exon():
    return Region(name="intronic", is_exon=False, pseudo=False)

@pytest.fixture
def mock_allele(functional_exon, pseudo_exon, non_exon):
    # Define a simple class to mimic the Allele behavior
    class MockAllele:
        @staticmethod
        def region(pos):
            if 10 <= pos <= 50:
                return functional_exon, pos
            elif 51 <= pos <= 100:
                return pseudo_exon, pos
            return non_exon, pos

    # Return an instance of the mock class
    return MockAllele()


def test_deletion_affects_functional_exon(mock_allele):
    assert identify_functional_deletion(15, "ACGT", mock_allele), "Deletion within a functional exon should return True"

def test_deletion_affects_only_pseudoexons(mock_allele):
    assert not identify_functional_deletion(60, "ACGT", mock_allele), "Deletion within only pseudoexons should return False"

def test_deletion_starts_in_functional_exon_ends_in_pseudoexon(mock_allele):
    assert identify_functional_deletion(45, "ACGT", mock_allele), "Deletion starting in functional exon and ending in pseudoexon should return True"

def test_deletion_starts_and_ends_in_same_non_functional_region(mock_allele):
    assert not identify_functional_deletion(110, "ACGT", mock_allele), "Deletion entirely within a non-functional region should return False"

# def test_deletion_spans_multiple_regions_including_functional_exon(mock_allele):
#     assert identify_functional_deletion(9, "ACGTACGT", mock_allele), "Deletion spanning multiple regions including a functional exon should return True"

# TODO
# def test_deletion_spans_multiple_regions_without_any_functional_exons(mock_allele):
#     assert not identify_functional_deletion(55, "ACGTACGT", mock_allele), "Deletion spanning multiple regions without any functional exons should return False"
