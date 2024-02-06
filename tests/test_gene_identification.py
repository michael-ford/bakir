import pytest
from kir_annotator.gene_identification import make_mapping_intervals, process_mapping, is_overlap_or_contained, update_interval_bounds, merge_overlapping_intervals

# MockAlignedSegment can either be defined here or imported if it's defined elsewhere
class MockAlignedSegment:
    def __init__(self, query_name, reference_start, reference_end):
        self.query_name = query_name
        self.reference_start = reference_start
        self.reference_end = reference_end
        self.is_unmapped = False
    def __eq__(self, other):
        if not isinstance(other, MockAlignedSegment):
            return NotImplemented
        return (self.query_name == other.query_name and 
                self.reference_start == other.reference_start and 
                self.reference_end == other.reference_end)

    def __hash__(self):
        return hash((self.query_name, self.reference_start, self.reference_end))

    def __repr__(self):
        return f"MockAlignedSegment({self.query_name}, {self.reference_start}, {self.reference_end})"
    
    def get_tag(self, tag):
        if tag == "NM":
            return 0
        else:
            raise KeyError(f"Tag {tag} not found")


make_interval_set = lambda interval: (interval[0], interval[1], set(interval[2]))
def test_make_mapping_intervals():
    mappings = [
        MockAlignedSegment("mapping1", 100, 200),
        MockAlignedSegment("mapping2", 150, 250),
        MockAlignedSegment("mapping3", 300, 400),
        MockAlignedSegment("mapping4", 350, 450),
        MockAlignedSegment("mapping5", 500, 600),
    ]


    intervals = make_mapping_intervals(mappings)

    assert len(intervals) == 3
    assert make_interval_set(intervals[0]) == (100, 250, set([mappings[0], mappings[1]]))  # Merged interval from mapping1 and mapping2
    assert make_interval_set(intervals[1]) == (300, 450, set([mappings[2], mappings[3]]))  # Merged interval from mapping3, mapping4
    assert make_interval_set(intervals[2]) == (500, 600, set([mappings[4]]))  # Single interval from mapping5

    mappings = [
        MockAlignedSegment("mapping1", 100, 200),
        MockAlignedSegment("mapping2", 150, 250),
        MockAlignedSegment("mapping3", 300, 400),
        MockAlignedSegment("mapping4", 350, 450),
        MockAlignedSegment("mapping5", 100, 450),  # Results in 1 interval
    ]

    intervals = make_mapping_intervals(mappings)
    print(intervals)

    assert len(intervals) == 1
    assert make_interval_set(intervals[0]) == (100, 450, set([mappings[0], mappings[1], mappings[2], mappings[3], mappings[4]]))  # Merged interval from all mappings   


def test_merge_overlapping_intervals():
    # Test case where two intervals are separate and should not be merged
    intervals = [
        (100, 200, ["mapping1"]),
        (300, 400, ["mapping2"])
    ]
    intervals = merge_overlapping_intervals(intervals)
    assert intervals == [(100, 200, ["mapping1"]), (300, 400, ["mapping2"])]

    # Test case where two intervals should be merged
    intervals = [
        (100, 200, ["mapping1"]),
        (150, 250, ["mapping2"]),
        (300, 400, ["mapping3"]),
        (350, 450, ["mapping4"]),
        (100, 450, ["mapping5"])  # Overlaps and connects the first two intervals
    ]
    intervals = merge_overlapping_intervals(intervals)
    print(intervals)
    assert len(intervals) == 1
    assert make_interval_set(intervals[0]) == (100, 450, set(["mapping1", "mapping2", "mapping5", "mapping3", "mapping4"]))


@pytest.mark.parametrize("start, end, mapping, expected", [
    (100, 200, MockAlignedSegment("mapping1", 150, 180), True),  # Mapping within interval
    (100, 200, MockAlignedSegment("mapping2", 50, 150), True),   # Mapping overlaps start
    (100, 200, MockAlignedSegment("mapping3", 150, 250), True),  # Mapping overlaps end
    (100, 200, MockAlignedSegment("mapping4", 50, 250), True),   # Mapping encompasses interval
    (100, 200, MockAlignedSegment("mapping5", 201, 300), False), # Mapping outside interval (after)
    (100, 200, MockAlignedSegment("mapping6", 50, 99), False)    # Mapping outside interval (before)
])
def test_is_overlap_or_contained(start, end, mapping, expected):
    assert is_overlap_or_contained(start, end, mapping) == expected

@pytest.mark.parametrize("start, end, mapping, expected", [
    (100, 200, MockAlignedSegment("mapping1", 150, 180), (100, 200)),  # Mapping entirely within interval
    (100, 200, MockAlignedSegment("mapping2", 50, 150), (50, 200)),    # Mapping expands start
    (100, 200, MockAlignedSegment("mapping3", 150, 250), (100, 250)),  # Mapping expands end
    (100, 200, MockAlignedSegment("mapping4", 50, 250), (50, 250))     # Mapping expands both ends
])
def test_update_interval_bounds(start, end, mapping, expected):
    assert update_interval_bounds(start, end, mapping) == expected

# Add more test functions as needed
