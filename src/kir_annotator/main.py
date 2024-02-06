from .io_utils import parse_arguments, read_sequence, read_database, write_output, extract_sample_info
from .mapping_utils import map_sequences
from .gene_identification import identify_genes
from .alignment_variants import identify_closest_allele, call_variants, align_to_wildtype, identify_functional_variants, adjust_sequence_from_variants
from .common import get_data_file_path, set_default_temp_dir
import logging
from collections import Counter, OrderedDict
import os
import mappy

logger = logging.getLogger("kir-annotator")

DATABASE_FASTA = 'kir_db.fasta'

def main(assembly_sequence_path: str, output_path: str = None, database_path: str = get_data_file_path('kir.pickle'), mapping_cache_path: str = None, temp_dir: str = None):
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
    database = read_database(database_path)

    mappings = map_sequences(assembly_sequence_path, get_data_file_path(DATABASE_FASTA), save_mapping_path=mapping_cache_path)
    gene_sequences = identify_genes(mappings, assembly_sequence, database)

    annotations = []
    for start, end, gene, gene_counts, assembly_seq, is_reverse in gene_sequences:
        if (sample_name, haplotype, gene, start) in to_fix:
            end = to_fix[(sample_name, haplotype, gene, start)]
            assembly_seq = str(assembly_sequence[mappings[0].reference_name][start:end]) if not is_reverse else mappy.revcomp(str(assembly_sequence[mappings[0].reference_name][start:end]))
        
        wildtype_sequence = database[gene].wildtype.seq
        
        assembly_seq, cigar_list = align_to_wildtype(assembly_seq, wildtype_sequence, seq_is_reversed=is_reverse)
        variants = call_variants(assembly_seq, cigar_list, wildtype_sequence)

        adjust_start, adjust_end = adjust_sequence_from_variants(variants, database[gene])

        if adjust_start or adjust_end:
            print(f"Testing adjustment of gene sequence for gene {gene} at interval ({start}, {end}). Adjustments: {adjust_start}, {adjust_end}. Is reverse: {is_reverse}")

            if not is_reverse:
                new_start = start+adjust_start
                new_end = end+adjust_end
                new_assembly_seq = str(assembly_sequence[mappings[0].reference_name][new_start:new_end])
            else:
                new_start = start-adjust_end
                new_end = end-adjust_start
                new_assembly_seq = mappy.revcomp(str(assembly_sequence[mappings[0].reference_name][new_start:new_end]))

            new_assembly_seq, new_cigar_list = align_to_wildtype(new_assembly_seq, wildtype_sequence, is_reverse)
            new_variants = call_variants(new_assembly_seq, new_cigar_list, wildtype_sequence)
            

            def size_variants(variants):
                return sum([len(op[3:]) if '>' not in op else 1 for _, op in variants])

            if size_variants(new_variants) < size_variants(variants):
                assembly_seq, cigar_list, is_reverse, variants, start, end = new_assembly_seq, new_cigar_list, new_is_reverse, new_variants, new_start, new_end
                print(f"Adjusted gene sequence for gene {gene} at interval ({start}, {end}). New sequence length: {len(assembly_seq)}")

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
            ('sequence', assembly_seq),
        ]), closest_alleles))


    write_output(annotations, output_path)


to_fix = {('HG00741', 'paternal', 'KIR3DP1', 373978): 376362}

def run() -> None:
    args = parse_arguments()
    main(assembly_sequence_path=args.sequence, mapping_cache_path=args.mapping_cache, output_path=args.output, temp_dir=args.temp_dir)
