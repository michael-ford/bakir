
import pysam
import logging
logger = logging.getLogger("kir-annotator")
from typing import List, Tuple
from .common import Seq
from .mapper_wrappers import Minimap2Wrapper
from collections import Counter

def identify_genes(mappings: List[pysam.AlignedSegment], assembly_sequence, database) -> List[Tuple[int, int, List[pysam.AlignedSegment]]]:
    """
    Identifies gene regions from mappings and groups overlapping mappings into intervals.

    Args:
        mappings (List[pysam.AlignedSegment]): Mappings from the assembly to the allele sequences.

    Returns:
        List[Tuple[int, int, List[pysam.AlignedSegment]]]: A list of intervals, each interval is a tuple with start, end, and a list of mappings within that interval.
    """
    intervals = make_mapping_intervals(mappings)
    
    updated_intervals = []
    for start, end, mappings in intervals:
        genes = [database[m.query_name.split('*')[0]].name for m in mappings]
        gene_counts = Counter(genes).most_common()

        if len(set(gene_counts)) > 1:
            logger.warning(f"More than 1 gene mapping for interval ({start}, {end}): {str(gene_counts)}")

        gene = genes[0]

        # remap to interval using wildtype to fine tune the sequence location
        wildtype_sequence = database[gene].wildtype.seq
        assembly_seq = str(assembly_sequence[mappings[0].reference_name][start:end])

        mapper = Minimap2Wrapper(params=f'-x map-ont')
        mapping_results = list(mapper.map([Seq('wildtype', wildtype_sequence)], Seq('assembly', assembly_seq)))

        start, end = start+mapping_results[0].reference_start, start+mapping_results[0].reference_end
        
        assembly_seq = str(assembly_sequence[mappings[0].reference_name][start:end])

        updated_intervals.append((start, end, gene, gene_counts, assembly_seq, mapping_results[0].is_reverse))
    
    return updated_intervals



def make_mapping_intervals(mappings: List[pysam.AlignedSegment]) -> List[Tuple[int, int, List[pysam.AlignedSegment]]]:
    """
    Groups mappings by unique mapping locations/intervals, merging overlapping mappings.

    Args:
        mappings (List[pysam.AlignedSegment]): List of pysam.AlignedSegment objects.

    Returns:
        List[Tuple[int, int, List[pysam.AlignedSegment]]]: List of intervals sorted by start position, each interval is a tuple of start, end, and mappings within the interval.
    """
    intervals = []  # Format: [[start, end, [mappings_within_interval, ...]]]

    for mapping in mappings:
        process_mapping(mapping, intervals)

    intervals = merge_overlapping_intervals(intervals)
    intervals.sort(key=lambda x: x[0])

    for i in range(len(intervals)):
        start, end, mappings = intervals[i]
        intervals[i] = (start, end, sorted(mappings, key=lambda x: x.get_tag('NM')))

    return intervals

def merge_overlapping_intervals(intervals: List[Tuple[int, int, List[pysam.AlignedSegment]]]) -> None:
    """
    Merges intervals that are overlapping or connected.

    Args:
        intervals (List[Tuple[int, int, List[pysam.AlignedSegment]]]): The list of intervals to merge.
    """
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    i = 0
    while i < len(intervals) - 1:
        current_interval = intervals[i]
        next_interval = intervals[i + 1]

        # Check if current and next intervals overlap or are connected
        if current_interval[1] >= next_interval[0]:
            # Merge the two intervals
            merged_interval = (
                current_interval[0], 
                max(current_interval[1], next_interval[1]), 
                current_interval[2] + next_interval[2]
            )

            # Update the list of intervals
            intervals[i] = merged_interval
            del intervals[i + 1]
        else:
            i += 1
    return intervals


def process_mapping(mapping: pysam.AlignedSegment, intervals: List[Tuple[int, int, List[pysam.AlignedSegment]]]) -> None:
    """
    Processes a single mapping and assigns it to the correct interval in the intervals list.

    Args:
        mapping (pysam.AlignedSegment): The mapping to process.
        intervals (List[Tuple[int, int, List[pysam.AlignedSegment]]]): The current list of intervals.
    """
    mapping_is_assigned = False

    for index, (start, end, interval_mappings) in enumerate(intervals):
        if is_overlap_or_contained(start, end, mapping):
            start, end = update_interval_bounds(start, end, mapping)
            interval_mappings.append(mapping)
            intervals[index] = (start, end, interval_mappings)
            mapping_is_assigned = True
            break  # Stop processing as mapping is assigned to an interval

    if not mapping_is_assigned:
        intervals.append((mapping.reference_start, mapping.reference_end, [mapping]))


def is_overlap_or_contained(start: int, end: int, mapping: pysam.AlignedSegment) -> bool:
    """
    Checks if a mapping overlaps with or is contained within an interval.

    Args:
        start (int): Start of the interval.
        end (int): End of the interval.
        mapping (pysam.AlignedSegment): The mapping to check.

    Returns:
        bool: True if the mapping overlaps with or is contained within the interval, False otherwise.
    """
    return (start <= mapping.reference_start <= end) or (start <= mapping.reference_end <= end) or \
           (mapping.reference_start < start and mapping.reference_end > end)


def update_interval_bounds(start: int, end: int, mapping: pysam.AlignedSegment) -> Tuple[int, int]:
    """
    Updates the bounds of an interval based on the mapping position.

    Args:
        start (int): Current start of the interval.
        end (int): Current end of the interval.
        mapping (pysam.AlignedSegment): The mapping to use for updating the interval bounds.

    Returns:
        Tuple[int, int]: Updated interval bounds.
    """
    new_start = min(start, mapping.reference_start)
    new_end = max(end, mapping.reference_end)
    return new_start, new_end


