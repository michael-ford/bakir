import re
import parasail
import mappy
from typing import Tuple, List, Dict
from .common import Allele, Gene
from Bio.Seq import Seq
from Bio import BiopythonWarning
import warnings
import logging
logger = logging.getLogger("kir-annotator")

def get_cigar(alignment):
    cigar_list = re.split(r'([A-Z=])', alignment.cigar.decode.decode())[:-1]
    return [(int(cigar_list[i]), cigar_list[i+1]) for i in range(0, len(cigar_list), 2)]


def align_to_wildtype(assembly_seq: str, wildtype_seq: str, seq_is_reversed: bool) -> Tuple[str, List[Tuple[int, str]]]:
    """
    Aligns a gene sequence from the assembly to the wildtype sequence of the closest allele.

    Args:
        assembly_seq (str): The gene sequence extracted from the assembly.
        wildtype_seq (str): The wildtype sequence of the closest allele.

    Returns:
        Tuple[str, List[Tuple[int, str]]]: The assembly seq as it was aligned (reverse compl as needed), and a list tuples of (pos, op) representing the CIGAR string. op is defined as X is a mismatch, I is an insertion, D is a deletion and = is a match.
    """
    if not assembly_seq or not wildtype_seq:
        raise ValueError("Empty sequence")

    def perform_alignment(seq1, seq2):
        parasail_matrix = parasail.matrix_create("ACGT", 2, -8)
        gap_open_penalty = 12
        gap_extend_penalty = 2
        return parasail.sg_qx_trace_scan_32(seq1, seq2, gap_open_penalty, gap_extend_penalty, parasail_matrix)

    if seq_is_reversed:
        assembly_seq = mappy.revcomp(assembly_seq)
    alignment = perform_alignment(wildtype_seq, assembly_seq)
    cigar_list = get_cigar(alignment)
    
    return assembly_seq, cigar_list


def call_variants(aligned_seq: str, cigar_list: List[Tuple[int, str]], wildtype_seq: str) -> set:
    """
    Calls variants from an alignment.

    Args:
        aligned_seq (str): The aligned sequence.
        wildtype_seq (str): The wildtype sequence for comparison.
        cigar_list (List[Tuple[int, str]]): list tuples of (pos, op) representing the CIGAR string. op is defined as X is a mismatch, I is an insertion, D is a deletion and = is a match.

    Returns:
        set: A set of variants represented as tuples (position, reference, alternate).
    """
    variants = set()
    wildtype_index, aligned_index = 0, 0

    for (size, operation) in cigar_list:
        if operation == "D":  # Deletion
            variants.add((wildtype_index, "ins" + aligned_seq[aligned_index:aligned_index + size]))
            aligned_index += size
        elif operation == "I":  # Insertion
            variants.add((wildtype_index, "del" + wildtype_seq[wildtype_index:wildtype_index + size]))
            wildtype_index += size
        else:
            if operation == "X":  # Mismatch
                for i in range(size):
                    variants.add((wildtype_index + i, f"{wildtype_seq[wildtype_index + i]}>{aligned_seq[aligned_index + i]}"))
            aligned_index += size
            wildtype_index += size

    return variants

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



def identify_functional_variants(variants: List[Tuple[int, str]], wildtype: Allele) -> List[Tuple[int, str]]:
    """
    Identifies functional variants among a list of variants.

    Args:
        variants (List[Tuple[int, str]]): List of variants represented as tuples (position, variant_str). variant_str is defined as ins{sequence} for insertions, del{sequence} for deletions and {wildtype}>{assembly} for mismatches.
        wildtype (Allele): The wildtype allele for reference.

    Returns:
        List[Tuple[int, str]]: A list of functional variants.
    """
    func_variants = []
    for pos, op in variants:
        region, relative_pos = wildtype.region(pos)
        if region and region.is_exon and not region.pseudo:
            if op.startswith("ins") or op.startswith("del"):
                func_variants.append((pos, op))
            else:
                assert wildtype.seq[pos] == op[0]
                seq = list(region.seq)
                seq[relative_pos] = op[2]  # apply SNP
                prot = "".join(
                    i.seq if i.name != region.name else "".join(seq)
                    for i in wildtype.regions.values()
                    if i.is_exon and not i.pseudo
                )
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter("always")
                    prot = str(Seq(prot).translate())
                if w:
                    for warning in w:
                        if isinstance(warning.message, BiopythonWarning):
                            logging.warning(f"Error translating sequence for variant {pos} {op} on wildtype {wildtype.gene.name+'*'+wildtype.name}: {str(warning.message)}")
                if prot != wildtype.protein:
                    func_variants.append((pos, op))
    return func_variants


def identify_closest_allele(variants: List[Tuple[int, str]], functional_variants: List[Tuple[int, str]], gene: Gene) -> Dict:
    """
    Identifies the closest matching allele based on shared functional variants.

    Args:
        variants (List[Tuple[int, str]]): List of variants found in the assembly gene sequence.
        gene (Gene): The gene object containing allele information.

    Returns:
        Dict: Information about allele similarity and differences.
    """

    def jaccard_distance(set1: set, set2: set) -> float:
        if set1 == set2:
            return 0
        return 1 - (len(set1 & set2) / len(set1 | set2)) if set1 | set2 else 1

    def get_allele_functional_variants(allele: Allele) -> set:
        return {(pos, op) for pos, op in allele.mutations if (pos, op) in allele.gene.functional}

    functional_variants = set(functional_variants)
    variants = set(variants)
    
    allele_similarity = []

    for allele in gene.alleles.values():
        allele_func_variants = get_allele_functional_variants(allele)
        similarity_metrics = {
            "closest allele": allele.gene.name+'*'+allele.name,
            "common mut": allele.mutations & set(variants),
            "new mut": set(variants) - allele.mutations,
            "missing mut": allele.mutations - set(variants),
            "missing functional mut": allele_func_variants - functional_variants,
            "new functional mut": functional_variants - allele_func_variants,
            "jaccard functional distance": jaccard_distance(allele_func_variants, functional_variants),
            "jaccard distance": jaccard_distance(allele.mutations, set(variants)),
        }
        allele_similarity.append(similarity_metrics)

    return sorted(allele_similarity, key=lambda x: (x["jaccard functional distance"], x['jaccard distance']))

