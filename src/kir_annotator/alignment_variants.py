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

def perform_alignment(seq1, seq2):
    parasail_matrix = parasail.matrix_create("ACGT", 2, -8)
    gap_open_penalty = 12
    gap_extend_penalty = 2
    return parasail.sg_qx_trace_scan_32(seq1, seq2, gap_open_penalty, gap_extend_penalty, parasail_matrix)


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
        if op.startswith("del"):
            func_variants.extend(identify_exon_deletions(pos, op, wildtype))
        elif region and region.is_exon:
            if op.startswith("ins"):
                func_variants.append((pos, op))
            else:
                assert wildtype.seq[pos] == op[0], f"Reference mismatch at position {pos} for variant {op} on wildtype {wildtype.gene.name+'*'+wildtype.name}"

                if wildtype.gene.name == 'KIR3DP1' and pos == 2114 and op == 'C>T':
                    func_variants.append((pos, op))
                    continue

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

def identify_exon_deletions(pos: int, op: str, wildtype: Allele) -> List[Tuple[int, str]]:
    """
    Identifies exon deletions based on a mutation operation.

    This function calculates the actual deletion impact on exons considering the operation
    start position and length. It handles full and partial deletions within exons while considering
    whether the exon is functional (non-pseudo).

    Args:
        pos (int): The starting position of the deletion.
        op (str): The deletion operation string, indicating the sequence being deleted.
        wildtype (Allele): The wildtype allele object with exon regions information.

    Returns:
        List[Tuple[int, str]]: A list of tuples where each tuple represents the start position
                               and the deletion string for each affected exon.
    """
    end_pos = pos + len(op[3:]) - 1  # Adjust based on the actual format of 'op'

    custom_deletions = []
    if wildtype.gene.name in ['KIR3DP1', 'KIR3DL3']:
        custom_deletions = custom_deletions_fix(pos, op, wildtype)
        if custom_deletions:
            return custom_deletions

    # Get all affected regions by the deletion
    affected_regions = identify_deletion_regions(pos, op, wildtype)

    op_seq = op[3:]
    exon_deletions = []  # Store deletions that occur within exons
    for region in affected_regions:
        if region.is_exon:
            if len(region.seq) == len(op_seq) and region.seq != op_seq:
                offset = pos - region.start
                if offset < 0:
                    if region.seq == op_seq[-offset:]+op_seq[:-offset]:
                        exon_deletions.append((region.start, f"del{region.seq}"))
                        break
                elif offset > 0:
                    if region.seq == op_seq[-offset:]+op_seq[:-offset]:
                        exon_deletions.append((region.start, f"del{region.seq}"))
                        break

            # Calculate overlap of the deletion with the exon
            overlap_start = max(pos, region.start)
            overlap_end = min(end_pos, region.end)
            overlap_seq = region.seq[overlap_start - region.start:overlap_end - region.start + 1]

            if overlap_seq:  # Check if there is an overlapping sequence
                exon_deletions.append((overlap_start, f"del{overlap_seq}"))

    return exon_deletions + custom_deletions

def custom_deletions_fix(pos: int, op: str, wildtype: Allele) -> List[Tuple[int, str]]:
    custom_deletions_3DP1 = {0: 
                        (335,
                        'delGGGGATGGAGATCTGGGCCCAGAGGTGGAGATATAGGCCTGGAGGTGGAGTTATGGGCCTGGAGTGGAGATCTGGGCCTGGAGTGGATATATGGGCCTGGAGATGGAGTGATGGGCCTAGAAGTGGAGATCTGGGTCTGGAGTGGAGATATGGGCCTGGAGGTGGAGATATGGGCCTGGAGTGGAGATCTGGGCCTGGAGTGGAGATAGGAACCTGGAGGGGAGATATGAGCCTGGAGTGAAGATATTGGCCTGGGATGGAGATATGGGCCTGGAGTGGAGACATGGGCCTGGAGGTGGAGATATGGGCCTGGAGGTGGAGACATGGGCCTAGAGGTGGATATCTGGGCCTGGAGTGGACATATGGGCCTAGGATGGAGATATGGGCCTGGGTGTGGAGATATGGGCTTGGGGTGGAGATATGGGCCTGGATTGGAGATATGGGTCTAGGGTGGAAATATTGGCCTGGAGTGGAGATATGGGCCTGGAGTGGAGATATGGGCTTGGGGTGGGGATAGGGGCCTGGGGTGCGGATATGGGCCTGCAGGCTGGGTCTCTACACAGCCGACAGCCCTGTTCTTGGGTGCAGGCTGGCACTGAGGGTGAGTTTCCCTTCAGCCCAGCAAGGGCCTGGCTACCAAGACTCACAGCCCAGTGGGGGCAGCAAGGGAGTCCTGGTTTGCCTGCAGATGGATGGTCCATCATGATCTTTCTTTCCAGGGTTCTTCTTGCTGCAGGGGGCCTGGACACATGAGGGTGAGTCCTTCTCCAAACCTTCGGGTGTCATCTCCCCACATAAGAGGATTTTCCTGAAACAGGAGGGAAGCCCGGTGGGGGATTTTCTTATAAACAAGGATGAGGAGACCCTGGGGTGCTCAGCCCACAGTTCCGACCTTGCCCTCCCCAGCCTTCCTTTCCCTTGGCTGAGTCAGGTTCTGTGGGAACCCGGGAGGGTAGACTGGGGTCCTCCAAGCTGGGCTGTGCGGCTGGGATGTGGTGTCACTGGCAGAGGAAGGGAGCAAAGCAGTGCTAGGAACAGCAGGCCTCTGAGGACAAAGGTGTAACTCACACCCTCCAGCGTTTCCATGACGGTAGGGGCTGCAGTGTGGCTGCTGTCATTCTACCTCAGAGGTGGGGGAACCCCAGCCAGGGCCCTGACCTTCCAAATCCTCTGTTGGGGGCTCAGTTGTGTATTGTGGTTCACACATTGGCTGATATTCCATTCACAAAGAACATGCCCTCGACTCCATGTCTATTTGTGTTGTTTTATGTGAGTAATCTTGCAGGATTAAAATCTAGTAGGAGTCCCTTACTCAGCACTTGCTCAAAGTTCTCAGCTGACACTTTTGTTGTAGAGAGACGCCAAGTCTATGCGGGGTGGGTCCTTCCTGTAGCCCTGGGCACCCAGGTGTGGTAGGAGCCTTAGAAAGTGGAAATGGGAGAATCTTCTGACACGTGGAGGGAGGGGCGGCTC')}
    custom_deletions_3DL3 = {(0,
  'delGTATGAGAGATTGGATCTGAGACGTGTTTTGAGTTGGTTATAGTGAAGGATGCAAGGTGTCAATTCTAGTTGGAACAATTTCCAGGAAGCCATGTTCTGCTCTTGACCAAACAGCCACTGGGCCTCATGCAAGGTAGAAATAGCCTGCATACGTCATCCTCCCATGATGTGGTCAGCATGTAAACTGCATGAGCCCCTCACAACATCCTGTGTGCTGCTGAACTGAGCTGGGGCGCAGCCGCCTGTCTGCACCGGCAGCACCATGTCGCTCATGGTCGTCAGCATGGCGTGTGTTGGTGAGTCCTGGAAGGGAATCGAGGGAGGGAGCGGTGGGGTGGAGATCTGGGCCTGGAGTGGAGATATGGGCCTGGAGTGGAGATATGGGCCTGGAGTGGAGATATAGGCCTGGAGTGGAGATATGGGCCTGGGGTGGAGATATGGGCCTGGAGTGGAGATATGGGCCTGGAACTGTAGATATGGGCCTGAAGTAGAGATATGGGCCTGGAGTAGAGATATGGGCCTGGAACTGTAGATATGGGCCTGGAGTGGAGATATTGGCTTGGAGTGCAGATATGGACCTGGAATTGAGATACGGGCCTGGAGGTGGAGATATGGGCCTAGAGTGGAGATATGGGCCTGGAGGTGGAGATATGGGCCTGGAACTGTAGATATGGGCCTGGAGTAGAGATATGGGCCTGGAGTGGAGATGTTGGCTTGGAGTGCAGATATGGGCCTGGAATGGAGACACGGGCCTGGAGGTGGAGATACAGGCCTGGAGGTGGAGATATGGGCCTGGAGTGTAGATATGGGCCTGGAGTAGAGATATAGGACAGAGGTGGAGATATAGGCCTGGAGTGGAGATATGGGCCTGGAGTAGAGATATAG'):
                        (262, 'insTCAT')}
    if wildtype.gene.name == 'KIR3DP1':
        if pos == 0 and len(op[3:]) >= 1473:
            return [custom_deletions_3DP1[pos]]
    if wildtype.gene.name == 'KIR3DL3':
        if (pos, op) in custom_deletions_3DL3:
            return [custom_deletions_3DL3[(pos, op)]]
    return []

def identify_deletion_regions(pos: int, op: str, wildtype: Allele) -> bool:
    """
    Determines if a deletion affects functional (non-pseudo) exon regions within an allele.

    Args:
        pos (int): The starting position of the deletion.
        op (str): The deletion operation string, indicating the sequence being deleted.
        wildtype (Allele): The wildtype allele object, which should have a method 'region' to get region info.

    Returns:
        bool: True if the deletion affects any functional exon regions, False otherwise.
    """
    start_region, _ = wildtype.region(pos)
    end_region, _ = wildtype.region(pos + len(op[3:]) - 1)

    # If the start and end regions are the same and it's not a functional exon, no impact
    if start_region.name == end_region.name:
        if start_region.is_exon:
            return [start_region]
        else:
            return []

    # Check if any region affected by the deletion is a functional exon
    affected_regions = []
    within_affected = False
    for region in wildtype.regions.values():
        # Start adding regions once the start_region is encountered
        if region.name == start_region.name:
            within_affected = True
        if within_affected:
            affected_regions.append(region)
        if region.name == end_region.name:
            break  # Stop once the end_region is passed

    return affected_regions

def identify_functional_deletion(pos: int, op: str, wildtype: Allele) -> bool:
    """
    Determines if a deletion affects functional (non-pseudo) exon regions within an allele.

    Args:
        pos (int): The starting position of the deletion.
        wildtype (Allele): The wildtype allele object, which should have a method 'region' to get region info.

    Returns:
        bool: True if the deletion affects any functional exon regions, False otherwise.
    """
    affected_regions = identify_deletion_regions(pos, op, wildtype)
    if any([region.is_exon and not region.pseudo for region in affected_regions]):
        return True
    return False



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
            "jaccard functional distance": normalized_jaccard_distance(allele_func_variants, gene, functional_variants),
            "jaccard distance": jaccard_distance(allele.mutations, set(variants)),
        }
        allele_similarity.append(similarity_metrics)

    return sorted(allele_similarity, key=lambda x: (x["jaccard functional distance"], x['jaccard distance']))

def normalized_jaccard_distance(allele_variants, gene, functional_variants: set) -> float:
    """
    Calculates the normalized Jaccard index between the called functional variants and the allele's variants,
    considering the gene's functional variants as the universe of possible variants.

    Args:
        allele: The allele object, expected to have a 'mutations' attribute.
        gene: The gene object, expected to have a 'functional' attribute which is a set of functional variant positions.
        functional_variants (set): A set of positions of functional variants called.

    Returns:
        float: The normalized Jaccard index.
    """
    # Extract called variant positions and allele variants within the gene's functional set
    called_var_pos = functional_variants

    # Calculate intersection: variants agree between allele and call, or both match the wildtype
    intersection = {var for var in gene.functional if (var in allele_variants) == (var in called_var_pos)}

    # Calculate union: all unique variants from gene's functional variants and called functional variants
    union = set(gene.functional.keys()) | called_var_pos

    # Return normalized Jaccard index: size of intersection divided by size of union
    return 1-(sum([calculate_variant_size(op) for p, op in intersection]) / sum([calculate_variant_size(op) for p, op in union])) if union else 1.0  # Return 1.0 if union is empty, meaning perfect match/no variation


def calculate_variant_size(op):
    """
    Calculates the size of a variant operation.

    Args:
        op (str): The variant operation string.

    Returns:
        int: The size of the variant operation.
    """
    if op.startswith("ins") or op.startswith("del"):
        return len(op[3:])
    else: # SNV
        return 1