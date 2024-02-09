from .io_utils import parse_arguments, read_sequence, read_database, write_output, extract_sample_info
from .mapping_utils import map_sequences
from .gene_identification import identify_genes
from .alignment_variants import identify_closest_allele, call_variants, align_to_wildtype, identify_functional_variants
from .common import get_data_file_path, set_default_temp_dir
from .annotation_utils import adjust_gene_sequence, apply_fixes_from_yaml
import logging
from collections import Counter, OrderedDict
import os
import mappy

logger = logging.getLogger("kir-annotator")

DATABASE_FASTA = 'kir_db.fasta'
DATABASE_PICKLE = 'kir_db.pickle'

def main(assembly_sequence_path: str, output_path: str = None, database_path: str = None, mapping_cache_path: str = None, temp_dir: str = None, yaml_fix_path: str = None):
    """
    The main driver function for the KIR annotator tool.

    This function orchestrates the process of gene sequence analysis, starting from reading the assembly sequence,
    mapping the sequences, identifying genes, and annotating the identified genes with relevant information such as
    functional variants and closest alleles.

    Args:
        assembly_sequence_path (str): The file path to the assembly sequence.
        output_path (str, optional): The file path for the output data. Defaults to None, which generates an output path based on the assembly_sequence_path.
        database_path (str, optional): The file path to the gene database. Defaults to kir.pickle in the data directory.
        mapping_cache_path (str, optional): The file path for caching mapping results. Defaults to None.
        temp_dir (str, optional): The file path to a temporary directory. Defaults to None.
        json_fix_path (str, optional): The file path to a JSON file containing manual fixes for gene sequences. Defaults to None.

    Returns:
        None: This function writes the output to a file and does not return any value.
    """

    if temp_dir:
        set_default_temp_dir(temp_dir)

    # Generate default output path if not provided
    if not output_path:
        output_path = os.path.splitext(os.path.basename(assembly_sequence_path))[0] + '_KIR-alleles'
    
    sample_name, haplotype = extract_sample_info(assembly_sequence_path)

    assembly_sequence = read_sequence(assembly_sequence_path)
    database = read_database(database_path if database_path else get_data_file_path(DATABASE_PICKLE))

    mappings = map_sequences(assembly_sequence_path, get_data_file_path(DATABASE_FASTA), save_mapping_path=mapping_cache_path)
    poss_gene_sequences = identify_genes(mappings, assembly_sequence, database)

    annotations = []
    for start, end, gene, gene_counts, poss_gene_seq, is_reverse, contig_sequence in poss_gene_sequences:
        
        start, end, assembly_seq = apply_fixes_from_yaml(sample_name, haplotype, gene, start, end, contig_sequence, yaml_fix_path)

        wildtype_sequence = database[gene].wildtype.seq
        
        poss_gene_seq, cigar_list = align_to_wildtype(poss_gene_seq, wildtype_sequence, seq_is_reversed=is_reverse)
        
        variants = call_variants(poss_gene_seq, cigar_list, wildtype_sequence)

        new_poss_gene_seq, new_cigar_list, new_is_reverse, new_variants, new_start, new_end = adjust_gene_sequence(start, end, database[gene], is_reverse, contig_sequence, wildtype_sequence, variants)
        poss_gene_seq = new_poss_gene_seq if new_poss_gene_seq else poss_gene_seq
        cigar_list = new_cigar_list if new_cigar_list else cigar_list
        is_reverse = new_is_reverse if new_is_reverse else is_reverse
        variants = new_variants if new_variants else variants
        start = new_start if new_start else start
        end = new_end if new_end else end

        functional_variants = identify_functional_variants(variants, database[gene].wildtype)
            
        closest_alleles = identify_closest_allele(variants, functional_variants, database[gene])

        annotations.append((OrderedDict([
            ('sample', sample_name),
            ('haplotype', haplotype),
            ('reference', mappings[0].reference_name),
            ('gene', gene),
            ('allele mapping counts', gene_counts),
            ('start', start),
            ('end', end),
            ('is reverse', is_reverse),
            ('sequence', poss_gene_seq),
        ]), closest_alleles))


    write_output(annotations, output_path)


to_fix = {('HG00741', 'paternal', 'KIR3DP1', 373978, 374305): {'end': 376362}}

def run() -> None:
    args = parse_arguments()
    main(assembly_sequence_path=args.sequence, 
        mapping_cache_path=args.mapping_cache, 
        output_path=args.output, 
        database_path=args.database,
        temp_dir=args.temp_dir,
        yaml_fix_path=args.fixes)
