# KIR Annotator 🧬

![Test Status](https://github.com/michael-ford/kir-annotator/actions/workflows/test.yml/badge.svg) [![Documentation Status](https://readthedocs.org/projects/your-project-name/badge/?version=latest)](https://your-project-name.readthedocs.io/en/latest/?badge=latest)

`kir-annotator` is a Python tool for gene sequence analysis, specifically focusing on the KIR (Killer-cell Immunoglobulin-like Receptor) genes. It automates the process of aligning sequences, identifying genes, and annotating them with relevant information such as functional variants and closest matching alleles.

## Features 🚀

- Sequence alignment and gene identification
- Variant calling and functional annotation
- Selection of best matching alleles through functional varaint identification
- Support for manual annotation fixes through YAML files
- Easy integration with CI/CD pipelines
- Comprehensive test suite with pytest

## Installation 📦

Clone the repository and navigate into the project directory:

```
git clone https://github.com/yourgithub/kir-annotator.git
cd kir-annotator
```
Install the environment with conda:
```
conda env create -f kir-annotator-env.yml
```

Install the package (preferably in a virtual environment):
```
pip install .
```

## Usage 🛠
To use kir-annotator, you need to provide the path to your assembly sequence file and the gene database file at a minimum. Additional options allow for customization of the output and other functionalities:

```
python -m src.kir_annotator.main sequence.fasta -d kir_db.fasta [-o output_file] [-m mapping_cache] [-t temp_dir] [-f fixes.yaml]
Command-Line Arguments
sequence: Path to the assembly sequence file.
-d, --database: Path to the gene database file.
-o, --output: (Optional) Path for the output file. Defaults to a name based on the sequence file name.
-m, --mapping_cache: (Optional) Path for the mapping cache file. If not provided, mapping results will not be cached.
-t, --temp_dir: (Optional) Path for the temporary directory. Uses the system default if not provided.
-f, --fixes: (Optional) Path to a YAML file containing manual annotation fixes.
```

## Development 🛠️
### Testing
We use pytest for testing. To run the tests:

`pytest tests/`

### CI/CD
This project is set up with CI/CD practices in mind. Automated tests are run on every push, ensuring code quality and functionality.

### Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue.
