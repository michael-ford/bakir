# BAKIR -  Biologically-informed Killer cell immunoglobulin-like receptor (KIR) gene annotation tool 🧬

![Test Status](https://github.com/michael-ford/bakir/actions/workflows/test.yml/badge.svg) [![Documentation Status](https://readthedocs.org/projects/bakir/badge/?version=latest)](https://bakir.readthedocs.io/en/latest/?badge=latest)

`bakir` is a Python tool for gene sequence analysis, specifically focusing on the KIR (Killer-cell Immunoglobulin-like Receptor) genes. It automates the process of aligning sequences, identifying genes, and annotating them with relevant information such as functional variants and closest matching alleles.

## Features 🚀

- Sequence alignment and gene identification
- Variant calling and functional annotation
- Selection of best matching alleles through functional varaint identification
- Support for manual annotation fixes through YAML files
- Easy integration with CI/CD pipelines
- Comprehensive test suite with pytest

## Installation 📦

### Native Installation
Clone the repository and navigate into the project directory:

```
git clone https://github.com/yourgithub/bakir.git
cd bakir
```
Install the environment with conda:
```
conda env create -f bakir-env.yml
conda activate bakir
```

Install the package (preferably in a virtual environment):
```
pip install .
```

### Using Docker
To run `bakir` using Docker, you first need to pull the Docker image from Docker Hub:

```
docker pull cdslsahinalp/kir-annotator:latest
```

After pulling the image, you can run `bakir` using:

``docker run -v $(pwd):/data cdslsahinalp/kir-annotator:latest [-o /data/OUTPUT] ... /data/sequence.fasta
```

Replace `$(pwd)` with the absolute path to the directory containing your data files if you are not running the command from within the data directory. This command mounts your current directory (`$(pwd)`) to `/data` inside the container.

### Using Singularity
For environments where Docker is not available, you can use Singularity to run the Docker image. First, pull the Docker image as a Singularity image:

```
singularity pull docker://cdslsahinalp/kir-annotator:latest
```

Then, run `bakir` using the Singularity image:

```
singularity exec --bind /path/to/data:/data kir-annotator_latest.sif bakir  [-o /data/OUTPUT] /data/sequence.fasta
```

As with Docker, replace `$(pwd)` with the absolute path to your data directory if necessary. This command binds your current directory to `/data` inside the Singularity container.

## Usage 🛠
To use `bakir`, you need to provide the path to your assembly sequence file and the gene database file at a minimum. Additional options allow for customization of the output and other functionalities:

```
usage: bakir [-h] [-d DATABASE] [-o OUTPUT] [-m MAPPING_CACHE] [-t TEMP_DIR] [-f FIXES] sequence

BAKIR: Biologically-informed Killer cell immunoglobulin-like receptor (KIR) gene annotation tool.

positional arguments:
  sequence              Path to the assembly sequence file.

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        Path to the gene database file.
  -o OUTPUT, --output OUTPUT
                        Optional: Path and prefix for the output file. If not provided, a default name based on the sequence file name will be used.
  -m MAPPING_CACHE, --mapping_cache MAPPING_CACHE
                        Optional: Path for the mapping cache file. If not provided, mapping results will not be cached.
  -t TEMP_DIR, --temp_dir TEMP_DIR
                        Optional: Path for the temporary directory. If not provided, the system default temporary directory will be used.
  -f FIXES, --fixes FIXES
                        Optional: Path to the YAML file containing manual annotation fixes. If not provided, no manual fixes will be applied.
```

## Development 🛠️
### Testing
We use pytest for testing. To run the tests:

`pytest tests/`

### CI/CD
This project is set up with CI/CD practices in mind. Automated tests are run on every push, ensuring code quality and functionality.

### Contributing
Contributions are welcome! Please feel free to submit a pull request or open an issue.
