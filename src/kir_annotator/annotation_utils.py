from .alignment_variants import align_to_wildtype, call_variants
from .io_utils import load_fixes
from typing import Tuple, List, Dict, Optional
from .common import Gene
import mappy

def apply_fixes_from_yaml(sample_name: str, haplotype: str, gene: str, start: int, end: int, contig_sequence, is_reverse: bool, yaml_path: Optional[str] = None) -> Tuple[int, int, str]:
    """
    Applies any necessary adjustments to the start and end positions and assembly sequence using fixes loaded from a YAML file.
    Adjustments are made based on a nested dictionary structure where each key maps to another dictionary specifying 'start' and 'end'.
    If 'start' or 'end' keys are not present, the respective variable values are not overwritten.

    Args:
        sample_name (str): The sample name.
        haplotype (str): The haplotype.
        gene (str): The name of the gene.
        start (int): The original start position of the gene sequence.
        end (int): The original end position of the gene sequence.
        contig_sequence: The sequence of the chromosome corresponding to start, end.
        is_reverse (bool): Indicates if the sequence is reversed.
        yaml_path (Optional[str]): The file path to the YAML file containing the fixes. Can be None.

    Returns:
        Tuple[int, int, str]: The possibly adjusted start and end positions and the (possibly adjusted) assembly sequence.
    """
    to_fix = {}
    if yaml_path is not None:
        to_fix = load_fixes(yaml_path)

    # Assuming the key format is converted appropriately when used with YAML
    key = (sample_name, haplotype, gene, start, end)  # Updated key to include 'end' for consistency with your data structure
    if key in to_fix:
        fix_details = to_fix.get(str(key), {})  # Convert tuple key to string if necessary
        start_adjustment = fix_details.get('start', start)
        end_adjustment = fix_details.get('end', end)

        adjusted_sequence = str(contig_sequence[start_adjustment:end_adjustment])
        if is_reverse:
            adjusted_sequence = mappy.revcomp(adjusted_sequence)  # Assuming mappy.revcomp is available

        return start_adjustment, end_adjustment, adjusted_sequence

    return start, end, str(contig_sequence[start:end])


def adjust_gene_sequence(start: int, 
                        end: int, 
                        gene: Gene, 
                        is_reverse: bool, 
                        contig_sequence: Dict[str, str], 
                        wildtype_sequence: str, 
                        variants: List[Tuple[str, str]]) -> Tuple[str, List[Tuple[str, str]], bool, List[Tuple[str, str]], int, int]:
    """
    Adjusts the gene sequence and recalculates variants based on the specified start and end adjustments.

    Args:
        start (int): The original start position of the gene sequence.
        end (int): The original end position of the gene sequence.
        gene (str): The name of the gene.
        is_reverse (bool): Indicates if the sequence is in reverse orientation.
        contig_sequence (Dict[str, str]): The dictionary containing assembly sequences.
        wildtype_sequence (str): The wildtype sequence of the gene.
        variants (List[Tuple[str, str]]): The list of variants identified in the original sequence.

    Returns:
        Tuple[str, List[Tuple[str, str]], bool, List[Tuple[str, str]], int, int]: The adjusted sequence, new variants list, updated is_reverse flag, original variants, and updated start and end positions.
    """
    adjust_start, adjust_end = adjust_sequence_from_variants(variants, gene)
    
    if adjust_start or adjust_end:
        print(f"Testing adjustment of gene sequence for gene {gene.name} at interval ({start}, {end}). Adjustments: {adjust_start}, {adjust_end}. Is reverse: {is_reverse}")

        new_start, new_end = calculate_new_positions(start, end, adjust_start, adjust_end, is_reverse)
        new_poss_gene_seq = retrieve_adjusted_sequence(contig_sequence, new_start, new_end, is_reverse)

        new_poss_gene_seq, new_cigar_list = align_to_wildtype(new_poss_gene_seq, wildtype_sequence, is_reverse)
        new_variants = call_variants(new_poss_gene_seq, new_cigar_list, wildtype_sequence)

        if size_variants(new_variants) < size_variants(variants):
            print(f"Adjusted gene sequence for gene {gene.name} at interval ({start}, {end}). New sequence length: {len(new_poss_gene_seq)}")
            return new_poss_gene_seq, new_cigar_list, is_reverse, new_variants, new_start, new_end

    return None, None, None, None, None, None


def adjust_sequence_from_variants(variants: List[Tuple[int, str]], gene: Gene) -> Tuple[int, int]:
    """
    Identifies insertions or deletions (indels) at the beginning or end of a gene sequence,
    shared across all alleles of the gene, and adjusts the sequence accordingly.

    Args:
        variants: A list of tuples, each containing a position (int) and an operation (str).
        gene: An instance of the Gene class, with attributes 'wildtype' and 'mutations'.

    Returns:
        A tuple of two integers representing the adjustments at the start and end of the sequence.

    Raises:
        ValueError: If multiple start or end indels are found in the variants.
    """
    def is_start_indel(pos: int, op: str) -> bool:
        """Check if the operation is a start indel."""
        return (op.startswith("ins") or op.startswith("del")) and pos == 0

    def is_end_indel(pos: int, op: str) -> bool:
        """Check if the operation is an end indel."""
        return (op.startswith("ins") or op.startswith("del")) and (pos == len(gene.wildtype.seq) or pos == len(gene.wildtype.seq) - len(op[3:]))

    start_adjust, end_adjust = 0, 0
    all_gene_mutation_types = [(p, op[:3]) for p, op in gene.mutations]

    for pos, op in variants:
        if is_start_indel(pos, op):
            if start_adjust:
                raise ValueError("Multiple start indels found in variants")
            if (pos, op) not in gene.mutations:
                adjust_value = len(op[3:])
                start_adjust = adjust_value if op.startswith("ins") else -adjust_value

        if is_end_indel(pos, op):
            if end_adjust:
                raise ValueError("Multiple end indels found in variants")
            if (pos, op) not in gene.mutations:
                adjust_value = len(op[3:])
                end_adjust = -adjust_value if op.startswith("ins") else adjust_value


    return start_adjust, end_adjust


def retrieve_adjusted_sequence(contig_sequence: str, new_start: int, new_end: int, is_reverse: bool) -> str:
    """
    Retrieves the adjusted sequence from the assembly.

    Args:
        contig_sequence (str-like): Sequence of contig where genne is present.
        new_start (int): The new start position.
        new_end (int): The new end position.
        is_reverse (bool): Indicates if the sequence is reversed.

    Returns:
        str: The adjusted gene sequence.
    """
    if not is_reverse:
        return str(contig_sequence[new_start:new_end])
    else:
        return mappy.revcomp(str(contig_sequence[new_start:new_end]))


def size_variants(variants: List[Tuple[str, str]]) -> int:
    """
    Calculates the size of the variants.

    Args:
        variants (List[Tuple[str, str]]): The list of variants.

    Returns:
        int: The calculated size of the variants.
    """
    return sum([len(op[3:]) if '>' not in op else 1 for _, op in variants])


def calculate_new_positions(start: int, end: int, adjust_start: int, adjust_end: int, is_reverse: bool) -> Tuple[int, int]:
    """
    Calculates new start and end positions after adjustments.

    Args:
        start (int): Original start position.
        end (int): Original end position.
        adjust_start (int): Adjustment to start position.
        adjust_end (int): Adjustment to end position.
        is_reverse (bool): Indicates if the sequence is reversed.

    Returns:
        Tuple[int, int]: The new start and end positions.
    """
    if not is_reverse:
        new_start = start + adjust_start
        new_end = end + adjust_end
    else:
        new_start = start - adjust_end
        new_end = end - adjust_start
    return new_start, new_end

