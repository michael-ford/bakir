from .mapper_wrappers import Minimap2Wrapper
import pysam, os
from typing import List

def map_sequences(assembly: str, database_fasta: str, save_mapping_path: str=None) -> List[pysam.AlignedSegment]:
    """
    Maps allele sequences to an assembly sequence using minimap2 and returns the alignments as pysam AlignedSegments.

    Args:
        assembly (str): The path to the assembly sequence file.
        database_fasta (str): The path to the FASTA file containing allele sequences.

    Returns:
        List[pysam.AlignedSegment]: A list of pysam AlignedSegments representing the mappings.
    """

    mapper = Minimap2Wrapper(params=f'-a -x map-ont')
    
    if save_mapping_path and os.path.exists(save_mapping_path):
        mapping_result = pysam.AlignmentFile(save_mapping_path, 'r')
    else:
        mapping_result = mapper.map(database_fasta, assembly, output_path=save_mapping_path)
    mappings = [x for x in mapping_result if not x.is_unmapped]
    
    return mappings
