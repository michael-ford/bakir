import pytest
from bakir.alignment_variants import identify_closest_allele


class MockAllele:
    def __init__(self, name, mutations):
        self.name = name
        self.mutations = mutations

class MockGene:
    def __init__(self, name, functional, alleles):
        self.name = name
        self.functional = dict([(x, '') for x in functional])
        self.alleles = dict([(x.name, x) for x in alleles])

# Define the functional mutations for the gene
functional_mutations = {(1, "A>G"), (2, "T>C"), (3, "G>A")}

# Create mock alleles
allele_no_functional = MockAllele("Allele1", set())
allele_two_functional = MockAllele("Allele2", set([(1, "A>G"), (2, "T>C"), (5, "insG")]))
allele_three_functional = MockAllele("Allele3", set([(1, "A>G"), (2, "T>C"), (3, "G>A"), (6, "G>T")]))
alleles = [allele_no_functional, allele_two_functional, allele_three_functional]

# Mock gene with three alleles
mock_gene = MockGene('Gene1', functional_mutations, alleles)
for a in alleles:
    a.gene = mock_gene

from bakir.alignment_variants import identify_closest_allele

def test_empty_variants():
    variants = []
    result = identify_closest_allele(variants, variants, mock_gene)
    # Expected: All alleles should have the same jaccard distance since there are no variants
    assert result[0]["closest allele"] == "Gene1*Allele1"

def test_variants_identical_to_allele2():
    variants = [(1, "A>G"), (5, "insG")]
    result = identify_closest_allele(variants, [variants[0]], mock_gene)
    print(result)
    assert result[0]["closest allele"] == "Gene1*Allele2"
    assert result[0]["jaccard distance"] < result[1]["jaccard distance"]
    assert result[1]["closest allele"] == "Gene1*Allele1"
    assert result[2]["closest allele"] == "Gene1*Allele3"

def test_variants_identical_to_allele3():
    variants = [(1, "A>G"), (3, "G>A"), (6, "G>T")]
    f_variants = [(1, "A>G"), (3, "G>A")]
    result = identify_closest_allele(variants, f_variants, mock_gene)
    for i in result:
        print(i)
    assert result[0]["closest allele"] == "Gene1*Allele3"
    assert result[0]["jaccard distance"] < result[1]["jaccard distance"]

def test_variants_wildtype_with_novel():
    variants = [(100, "A>G"), (200, "T>C"), (400, "C>A")]
    f_variants = [(100, "A>G"), (200, "T>C"), (5, "insG")]
    result = identify_closest_allele(variants, f_variants, mock_gene)
    for i in result:
        print(i)
    assert result[0]["closest allele"] == "Gene1*Allele1"
    assert result[0]["jaccard functional distance"] < result[1]["jaccard functional distance"]
    assert result[1]["closest allele"] == "Gene1*Allele2"
    assert result[2]["closest allele"] == "Gene1*Allele3"