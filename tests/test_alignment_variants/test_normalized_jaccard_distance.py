# import pytest
# from collections import namedtuple
# from bakir.alignment_variants import normalized_jaccard_distance

# # Create mock structures
# Allele = namedtuple('Allele', ['mutations'])
# Gene = namedtuple('Gene', ['functional'])

# # Define pytest fixtures for common setups
# @pytest.fixture
# def mock_gene():
#     return Gene(functional={10, 20, 30, 40, 50})

# @pytest.fixture
# def mock_allele_full_overlap(mock_gene):
#     return Allele(mutations=mock_gene.functional)

# @pytest.fixture
# def mock_allele_partial_overlap():
#     return Allele(mutations={20, 30})

# @pytest.fixture
# def mock_allele_no_overlap():
#     return Allele(mutations={60, 70})

# def test_with_no_overlap(mock_gene, mock_allele_no_overlap):
#     functional_variants = {80, 90}
#     result = normalized_jaccard_distance(mock_allele_no_overlap, mock_gene, functional_variants)
#     assert result == 1, "Expected Jaccard index to be 0 with no overlap"

# def test_with_perfect_overlap(mock_gene, mock_allele_full_overlap):
#     functional_variants = mock_gene.functional
#     result = normalized_jaccard_distance(mock_allele_full_overlap, mock_gene, functional_variants)
#     assert result == 0, "Expected Jaccard index to be 1 with perfect overlap"

# def test_with_partial_overlap(mock_gene, mock_allele_partial_overlap):
#     functional_variants = {20, 30, 60}
#     result = normalized_jaccard_distance(mock_allele_partial_overlap, mock_gene, functional_variants)
#     assert result == 1- (2/6), "Expected Jaccard index to reflect partial overlap"

# def test_with_allele_variants_only(mock_gene, mock_allele_partial_overlap):
#     functional_variants = {60, 70}  # Non-overlapping with gene functional variants
#     result = normalized_jaccard_distance(mock_allele_partial_overlap, mock_gene, functional_variants)
#     assert result == 1- (2 / 7), "Expected Jaccard index to reflect overlap with allele variants only"

# def test_with_called_variants_only(mock_gene):
#     allele = Allele(mutations={60, 70})  # Non-overlapping with gene functional variants
#     functional_variants = {30, 40}  # Overlapping with gene functional variants
#     result = normalized_jaccard_distance(allele, mock_gene, functional_variants)
#     assert result == 1 -(2 / 7), "Expected Jaccard index to reflect overlap with called variants only"

# def test_with_matching_wildtype(mock_gene):
#     allele = Allele(mutations={})  # Partially overlapping with gene functional variants
#     functional_variants = {90, 70}  # Partially overlapping with gene functional variants
#     result = normalized_jaccard_distance(allele, mock_gene, functional_variants)
#     assert result == 1- (5 / 7), "Expected Jaccard index to reflect partial matching including wildtype"

# def test_with_empty_sets():
#     gene = Gene(functional=set())
#     allele = Allele(mutations=set())
#     functional_variants = set()
#     result = normalized_jaccard_distance(allele, gene, functional_variants)
#     assert result == 1, "Expected Jaccard index to be 1 when all sets are empty"

# def test_with_no_functional_variants():
#     gene = Gene(functional=set())  # No functional variants
#     allele = Allele(mutations={10, 20, 30})
#     functional_variants = {40, 50, 60}
#     result = normalized_jaccard_distance(allele, gene, functional_variants)
#     assert result == 1, "Expected Jaccard index to be 1 when there are no gene functional variants"
