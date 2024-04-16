import pytest
from collections import OrderedDict, namedtuple
from bakir.annotation_utils import filter_annotations
import logging

# Define the namedtuple for the mock database
Gene = namedtuple('Gene', ['wildtype'])
Wildtype = namedtuple('wildtype', ['seq'])

# Fixture for mock annotations
@pytest.fixture
def mock_annotations():
    return [
        (OrderedDict([
            ('sample', 'sample1'), 
            ('gene', 'gene1'), 
            ('start', 100), 
            ('end', 300), 
            ('sequence', 'ACTG')
        ]), 'allele1'),
        (OrderedDict([
            ('sample', 'sample2'), 
            ('gene', 'gene2'), 
            ('start', 200), 
            ('end', 400), 
            ('sequence', 'GTCAG')
        ]), 'allele2'),
    ]

# Fixture for mock database
@pytest.fixture
def mock_database():
    return {
        'gene1': Gene(wildtype=Wildtype(seq='A'*50)),  # 50 As
        'gene2': Gene(wildtype=Wildtype(seq='C'*50)),  # 50 Cs
    }

def test_non_overlapping_annotations(mock_annotations, mock_database):
    filtered, incomplete = filter_annotations(mock_annotations, mock_database, min_cov=0.8)
    assert len(filtered) == len(mock_annotations), "All non-overlapping annotations should be retained"
    assert not incomplete, "There should be no incomplete gene copies when annotations are non-overlapping"

def test_one_annotation_contained_within_another(mock_annotations, mock_database):
    # Modify the second annotation to be contained within the first
    mock_annotations[1] = (mock_annotations[1][0].copy(), mock_annotations[1][1])
    mock_annotations[1][0]['start'], mock_annotations[1][0]['end'] = 150, 250  # Contained within the first
    filtered, incomplete = filter_annotations(mock_annotations, mock_database)
    assert len(filtered) == 1 and filtered[0][0]['sample'] == 'sample1', "Only the containing annotation should be retained"
    assert not incomplete, "There should be no incomplete gene copies when one annotation is contained within another"

def test_multiple_nested_annotations(mock_annotations, mock_database):
    # Add a third annotation that is contained within the second
    mock_annotations.append((OrderedDict([
        ('sample', 'sample3'), 
        ('gene', 'gene3'), 
        ('start', 220), 
        ('end', 240), 
        ('sequence', 'GGG')
    ]), 'allele3'))
    mock_database['gene3'] = Gene(wildtype=Wildtype(seq='G'*20))  # 50 Gs
    filtered, incomplete = filter_annotations(mock_annotations, mock_database)
    assert len(filtered) == len(mock_annotations) - 1, "Only the outermost annotations should be retained"
    assert not incomplete, "There should be no incomplete gene copies when handling nested annotations"

def test_identical_annotations(mock_annotations, mock_database):
    # Add an identical annotation
    mock_annotations.append((mock_annotations[0][0].copy(), mock_annotations[0][1]))
    filtered, incomplete = filter_annotations(mock_annotations, mock_database)
    assert len(filtered) == len(mock_annotations) - 1, "Only one of the identical annotations should be retained"
    assert not incomplete, "There should be no incomplete gene copies with identical annotations"

def test_all_annotations_contained(mock_annotations, mock_database):
    # Modify annotations so each is contained within the previous
    mock_annotations[1] = (mock_annotations[1][0].copy(), mock_annotations[1][1])
    mock_annotations[1][0]['start'], mock_annotations[1][0]['end'] = 150, 250
    mock_annotations.append((OrderedDict([
        ('sample', 'sample3'), 
        ('gene', 'gene3'), 
        ('start', 170), 
        ('end', 230), 
        ('sequence', 'GGGG')
    ]), 'allele3'))
    mock_database['gene3'] = Gene(wildtype=Wildtype(seq='G'*50))  # 50 Gs
    filtered, incomplete = filter_annotations(mock_annotations, mock_database)
    assert len(filtered) == 1, "Only the outermost annotation should be retained when all are contained"
    assert not incomplete, "There should be no incomplete gene copies when all annotations are contained"

def test_coverage_filter(mock_annotations, mock_database):
    # Modify an annotation to have insufficient coverage
    mock_annotations[1][0]['start'] = 0  # Adjust to change coverage
    mock_annotations[1][0]['end'] = 1  # Adjust to change coverage
    filtered, incomplete = filter_annotations(mock_annotations, mock_database)
    print(filtered)
    print(incomplete)
    assert len(incomplete) == 1 and incomplete[0]['sample'] == 'sample2', "Annotations below the coverage threshold should be marked incomplete"

def test_empty_annotations_list(mock_database):
    filtered, incomplete = filter_annotations([], mock_database)
    assert not filtered and not incomplete, "Should return empty lists when no annotations are provided"

def test_annotations_with_insufficient_coverage(mock_annotations, mock_database):
    # Set all annotations to have insufficient coverage
    for ann in mock_annotations:
        ann[0]['end'] = ann[0]['start'] + 10  # Dramatically reduce length
    filtered, incomplete = filter_annotations(mock_annotations, mock_database, min_cov=0.8)
    assert not filtered and len(incomplete) == len(mock_annotations), "All annotations should be marked incomplete due to low coverage"
