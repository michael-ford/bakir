import pytest
from kir_annotator.alignment_variants import call_variants


def test_call_variants_mismatch():
    assembly_seq = "AGCTT"
    wildtype_seq = "AGCCT"  # Mismatch at position 3
    cigar_list = [(3, "="), (1, "X"), (1, "=")]

    variants = call_variants(assembly_seq, cigar_list, wildtype_seq)
    assert variants == {(3, "C>T")}

def test_call_variants_insertion():
    assembly_seq = "AGCGTA"
    wildtype_seq = "AGCTA"  # Insertion at position 3
    cigar_list = [(3, "="), (1, "D"), (2, "=")]

    variants = call_variants(assembly_seq, cigar_list, wildtype_seq)
    assert variants == {(3, "insG")}

def test_call_variants_homo_insertion():
    assembly_seq = "AGCGGGTA"
    wildtype_seq = "AGCTA"  # Insertion at position 3
    cigar_list = [(3, "="), (3, "D"), (2, "=")]

    variants = call_variants(assembly_seq, cigar_list, wildtype_seq)
    assert variants == {(3, "insGGG")}


def test_call_variants_deletion():
    assembly_seq = "AGTA"
    wildtype_seq = "AGCTA"  # Deletion at position 3
    cigar_list = [(2, "="), (1, "I"), (2, "=")]

    variants = call_variants(assembly_seq, cigar_list, wildtype_seq)
    assert variants == {(2, "delC")}

def test_call_variants_insertion_start():
    assembly_seq = "ACGTAAGCA"
    wildtype_seq = "AAGCA"  
    cigar_list = [(4, 'D'), (5, '=')]

    variants = call_variants(assembly_seq, cigar_list, wildtype_seq)
    assert variants == {(0, "insACGT")}

def test_call_variants_deletion_start():
    assembly_seq = "AGCT"
    wildtype_seq = "TTTAGCT"  
    cigar_list = [(3, 'I'), (4, '=')]

    variants = call_variants(assembly_seq, cigar_list, wildtype_seq)
    assert variants == {(0, "delTTT")}

def test_call_variants_insertion_end():
    assembly_seq = "ACGTAGCG" 
    wildtype_seq = "ACGTA"
    cigar_list = [(5, '='), (3, 'D')]

    variants = call_variants(assembly_seq, cigar_list, wildtype_seq)
    assert variants == {(5, "insGCG")}

def test_call_variants_deletion_end():
    assembly_seq = "AGCT"
    wildtype_seq = "AGCTAAA"
    cigar_list = [(4, '='), (3, 'I')]

    variants = call_variants(assembly_seq, cigar_list, wildtype_seq)
    assert variants == {(4, "delAAA")}
