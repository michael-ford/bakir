import sys
import os
import pickle
from typing import Dict, Tuple, List, Optional
from .common import Gene
import argparse
import pyfastx
import yaml
from collections import OrderedDict
from yaml.loader import SafeLoader


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for the KIR annotator tool.

    This function sets up command-line arguments for the script, enabling users to specify the paths for the
    assembly sequence file, the gene database file, and optionally the output path, the mapping cache path,
    and the path to a YAML file containing manual annotation fixes.

    Returns:
        argparse: Argparse argument that needs .parse_args() call.
            - sequence (str): Path to the assembly sequence file.
            - database (str): Path to the gene database file.
            - output (str, optional): Path for the output file. Defaults to None.
            - mapping_cache (str, optional): Path for the mapping cache file. Defaults to None.
            - temp_dir (str, optional): Path for the temporary directory. Defaults to None.
            - fixes (str, optional): Path to the YAML file containing manual annotation fixes. Defaults to None.
    """
    parser = argparse.ArgumentParser(description="KIR Annotator Tool: A tool for gene sequence analysis.")

    parser.add_argument(
        'sequence',
        type=str,
        help='Path to the assembly sequence file.'
    )

    parser.add_argument(
        '-d', '--database',
        type=str,
        default=None,
        help='Path to the gene database file.'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Optional: Path and prefix for the output file. If not provided, a default name based on the sequence file name will be used.'
    )

    parser.add_argument(
        '-m', '--mapping_cache',
        type=str,
        default=None,
        help='Optional: Path for the mapping cache file. If not provided, mapping results will not be cached.'
    )

    parser.add_argument(
        '-t', '--temp_dir',
        type=str,
        default=None,
        help='Optional: Path for the temporary directory. If not provided, the system default temporary directory will be used.'
    )

    parser.add_argument(
        '-f', '--fixes',
        type=str,
        default=None,
        help='Optional: Path to the YAML file containing manual annotation fixes. If not provided, no manual fixes will be applied.'
    )

    return parser

def read_sequence(file: str) -> pyfastx.Fasta:
    """
    Reads a genomic sequence from a given file.
    
    Args:
        file (str): The path to the file containing the genomic sequence.

    Returns:
        str: The genomic sequence read from the file.
    """
    return pyfastx.Fasta(file)


def extract_sample_info(assembly_sequence_path: str) -> Tuple[str, str]:
    """
    Extracts sample name and haplotype from the assembly sequence file name.

    Args:
        assembly_sequence_path (str): The file path to the assembly sequence.

    Returns:
        Tuple[str, str]: The sample name and haplotype.
    """
    file_name = os.path.basename(assembly_sequence_path)
    sample_name, haplotype = file_name.split('.')[:2]
    return sample_name, haplotype



def load_pickle(file_path):
    """
    Loads a pickle file, ensuring that the 'kir_annotator' package is in the Python path.

    Args:
        file_path (str): Path to the pickle file.

    Returns:
        The unpickled object.
    """
    # Determine the directory of the kir_annotator package
    kir_annotator_dir = os.path.dirname(os.path.abspath(__file__))
    # Insert this directory into sys.path
    if kir_annotator_dir not in sys.path:
        sys.path.insert(0, kir_annotator_dir)

    # Now load the pickle file
    with open(file_path, 'rb') as file:
        obj = pickle.load(file)
    return obj


def read_database(file: str) -> Dict[str, Gene]:
    """
    Reads a database of genes and alleles from a given file.
    
    Args:
        file (str): The path to the file containing the database.

    Returns:
        Dict[str, Gene]: A dictionary with gene IDs as keys.
    """
    
    return load_pickle(file)[0]


def load_fixes(yaml_path: str) -> Dict:
    """
    Loads and caches the fixes dictionary from a YAML file.

    Args:
        yaml_path (str): The file path to the YAML file containing the fixes.

    Returns:
        Dict: The loaded dictionary of fixes.

    Function Attributes:
        cache (dict): A cache storing loaded dictionaries of fixes, keyed by their YAML file paths.
                      This attribute is dynamically attached to the function to retain state across calls.
    """
    if not hasattr(load_fixes, "cache"):
        load_fixes.cache = {}
    if yaml_path not in load_fixes.cache:
        with open(yaml_path, 'r') as file:
            load_fixes.cache[yaml_path] = yaml.load(file, Loader=SafeLoader)
    return load_fixes.cache[yaml_path]


def summarize_output_data(data: List[Dict[str, List[Dict[str, any]]]]) -> Tuple[List[Dict], str]:
    """
    Summarize output data by extracting and calculating relevant information.

    Args:
        data (List[Dict[str, List[Dict[str, any]]]]): A list of dictionaries, each containing annotation data.

    Returns:
        Tuple[List[Dict], str]: A tuple containing a list of summarized annotations and a formatted TSV string.
    """

    output_dicts = []
    output_tsv = ["sample\thaplotype\tgene\tallele\tstart\tend\tnew variants\tnew functional_variants\tsequence\n"]

    for annotation, closest_alleles in data:
        new_closest_alleles = []

        for allele in closest_alleles[:5]:
            new_closest_alleles.append(OrderedDict({
                "allele": allele.get("closest allele", ""),
                "new mut": len(allele.get("new mut", [])),
                "missing mut": len(allele.get("missing mut", [])),
                "new functional mut": allele.get("new functional mut", []),
                "missing functional mut": allele.get("missing functional mut", []),
                "jaccard functional distance": allele.get("jaccard functional distance", 0.0),
                "jaccard distance": allele.get("jaccard distance", 0.0),
            }))

        # Update annotation with new information
        annotation["closest allele"] = new_closest_alleles[0]["allele"] if new_closest_alleles else None
        annotation["new functional variants"] = closest_alleles[0]["new functional mut"] if closest_alleles else None
        annotation["missing functional variants"] = closest_alleles[0]["missing functional mut"] if closest_alleles else None
        annotation["new mut"] = closest_alleles[0]["new mut"] if closest_alleles else None
        annotation["missing mut"] = closest_alleles[0]["missing mut"] if closest_alleles else None
        annotation["top alleles"] = new_closest_alleles
        output_dicts.append(annotation)

        # Format data for TSV
        output_tsv.append(f"{annotation['sample']}\t{annotation['haplotype']}\t{annotation['gene']}\t{annotation['closest allele']}\t{annotation['start']}\t{annotation['end']}\t{len(closest_alleles[0].get('variants', []))}\t{len(closest_alleles[0].get('new functional variants', []))}\t{annotation['sequence']}\n")

    return output_dicts, ''.join(output_tsv)


def write_output(data: List[Dict], output_path: str) -> None:
    """
    Writes the analysis results to a file.
    
    Args:
        data (Dict[str, List[str]]): The data to write.
        file (str): The path to the output file.
    """

    yaml_data, tsv_data = summarize_output_data(data)

    with open(output_path+'.yaml', 'w') as f:
        yaml.dump(yaml_data, f)

    with open((output_path)+'.tsv', 'w') as f:
        f.write(tsv_data)
        
            
