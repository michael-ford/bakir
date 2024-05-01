import os
import sys
import time
from datetime import datetime
from dataclasses import dataclass
import os
from typing import Optional
import tempfile
import logging

# Define a format to include the file name and line number
log_format = '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s'

# Configure your logger
logging.basicConfig(level=logging.INFO, format=log_format)
logger = logging.getLogger("kir-annotator")


class Seq:
    def __init__(self, name, seq):
        self.id = name
        self.seq = seq

def get_data_file_path(filename: str) -> Optional[str]:
    """
    Constructs the path to a data file in the package's data directory.

    Args:
        filename (str): The name of the file in the data directory.

    Returns:
        Optional[str]: The full path to the data file, or None if the file does not exist.
    """
    # Get the directory of the current file (gene_identification.py)
    current_dir = os.path.dirname(__file__)

    # Construct the path to the data file
    data_file_path = os.path.join(current_dir, 'data', filename)

    # Check if the file exists
    if os.path.exists(data_file_path):
        return data_file_path
    else:
        logger.error(f"Error: The file {filename} does not exist.")
        return None


@dataclass
class Region:
    gene: str
    name: str
    start: int
    end: int
    seq: str
    partial: bool = False
    pseudo: bool = False

    @property
    def is_exon(self):
        return self.name.startswith("e")

    @property
    def is_intron(self):
        return self.name.startswith("i")

@dataclass
class Allele:
    gene: object
    name: str
    regions: dict
    ops: list  # list of differences from the wildtype allele
    protein: str

    enabled: bool
    func: set
    minor: set
    mutations: set
    keystones: set

    def __init__(self, gene, name, record):
        """Initialize the allele from an IMGT (kir.dat) record"""

        self.gene = gene
        self.name = name
        self.protein = None
        self.description = record.description
        self.regions = {}
        for f in record.features:
            if f.type == "CDS":
                if (protein := f.qualifiers.get("translation", None)):
                    self.protein = protein[0] + "X"
            elif f.type != "source":
                name = f.type.lower()
                if name in ["exon", "intron"]:
                    name = name[0]
                name += str(
                    f.qualifiers.get(
                        "number", ["_3" if f.type.lower() in self.regions else ""]
                    )[0]
                )
                assert name not in self.regions, (name, self.regions.keys())
                self.regions[name] = Region(
                    gene,
                    name,
                    int(f.location.start),
                    int(f.location.end),
                    str(record.seq[int(f.location.start) : int(f.location.end)]),
                    "partial" in f.qualifiers,
                    "pseudo" in f.qualifiers,
                )

        self.enabled = True
        self.func = set()
        self.minor = set()
        self.mutations = set()
        self.keystones = set()

    def parse_mutations(self):
        """Find all core mutations that modify the protein."""
        from Bio.Seq import Seq

        for pos, op in self.ops:
            self.gene.mutations.setdefault((pos, op), set()).add(self.name)
            if (pos, op) in self.gene.functional:
                continue
            r, rs = self.gene.wildtype.region(pos)
            if r and r.is_exon and not r.pseudo:  # ignore pseudo-exons!
                if op.startswith("ins") or op.startswith("del"):
                    # Exonic indels are always functionsl
                    self.gene.functional[pos, op] = op[:3]
                else:
                    # For SNPs, check if the protein is modified
                    seq = list(r.seq)
                    seq[rs] = op[2]  # apply the change
                    prot = "".join(  # translate to protein
                        i.seq if i.name != r.name else "".join(seq)
                        for i in self.gene.wildtype.regions.values()
                        if i.is_exon and not i.pseudo
                    )
                    prot = str(Seq(prot).translate())
                    if prot != self.gene.wildtype.protein:
                        pfx = os.path.commonprefix([prot, self.gene.wildtype.protein])
                        p = f"{self.gene.wildtype.protein[len(pfx) - 1]}{len(pfx)}{prot[len(pfx)]}"
                        self.gene.functional[pos, op] = p

    def region(self, i):
        """Return the region name and offset of gene location."""

        st = 0
        for r in self.regions.values():
            if st <= i < st + len(r.seq):
                return r, i - st
            st += len(r.seq)
        # assert False, i
        return None, None

    @property
    def seq(self):
        return "".join(s.seq for s in self.regions.values())

    @property
    def exons(self):
        return [r for r in self.regions.values() if r.is_exon]

    @property
    def introns(self):
        return [r for r in self.regions.values() if r.is_intron]

    @property
    def utrs(self):
        return [r for r in self.regions.values() if r.name.startswith("utr")]


class Gene:
    def __init__(self, gene, wildtype):
        self.name = gene
        self._wildtype = wildtype  # Wildtype allele
        self.alleles = {}
        self.functional = {}  # Functional mutations
        self.mutations = {}

        # The start position in the final chr19kir chromosome (generated later)
        self.ref_start = 0

    def add(self, allele, record):
        self.alleles[allele] = Allele(self, allele, record)

    @property
    def wildtype(self):
        return self[self._wildtype]

    def __getitem__(self, i):
        return self.alleles[i]

    @property
    def seq(self):
        return self.wildtype.seq

    def yaml(self):
        """Generate Aldy database YAML"""

        maps, exons = {}, []
        start = 0
        en = 1
        for r in self.wildtype.regions.values():
            name = r.name  # handle skip exons
            if name.startswith('e'):
                name = f"e{en}"
            if name.startswith('i'):
                name = f"i{en}"; en += 1
            maps[name] = [
                self.ref_start + start + 1,
                self.ref_start + start + len(r.seq) + 1,
            ]
            if r.is_exon and not r.pseudo:
                exons.append([start + 1, start + len(r.seq) + 1])
            start += len(r.seq)
        return {
            "name": self.name,
            "version": "fresh",
            "generated": str(datetime.now()),
            "alleles": {
                f"{self.name}*{a.name}": {
                    "mutations": [
                        [pos + 1 - (1 if op.startswith("ins") else 0), op, "-"]
                        + ([f] if (f := self.functional.get((pos, op))) else [])
                        for (pos, op) in sorted(a.ops)
                    ]
                }
                for a in self.alleles.values()
            },
            "structure": {
                "genes": [self.name],
                "regions": {"hg38": maps},
                "cn_regions": list(maps.keys()),
            },
            "reference": {
                "name": self.name,
                "mappings": {
                    "hg38": [
                        "19kir",
                        self.ref_start + 1,
                        self.ref_start + (l := len(self.seq)) + 1,
                        "+",
                        f"M{l}",
                    ]
                },
                "exons": exons,
                "seq": self.seq,
            },
        }


def bisect_left(a, x, lo=0, hi=None, *, key=None):
    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = len(a)

    if key is None:
        while lo < hi:
            mid = (lo + hi) // 2
            if a[mid] < x: lo = mid + 1
            else: hi = mid
    else:
        while lo < hi:
            mid = (lo + hi) // 2
            if key(a[mid]) < x: lo = mid + 1
            else:
                hi = mid
    return lo

class TimeInterval:
    def __init__(self, msg=""):
        self.start = time.time()
        self.msg = msg

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, *args):
        print(self.report(self.msg), file=sys.stderr)

    def report(self, msg="") -> str:
        import psutil, resource
        process = psutil.Process()
        mi = process.memory_full_info()
        mp = mi.rss / (1024 ** 2)
        ru = resource.getrusage(resource.RUSAGE_SELF)
        mx = ru.ru_maxrss / 1024
        msg = 'Block' if not self.msg else self.msg
        msg = f"[time] {msg} took {self.elapsed():.2f}s ({mp:,} MB; {mx:,} MB)"
        return msg

    def elapsed(self) -> float:
        return time.time() - self.start

def timing(msg: str = "") -> TimeInterval:
    return TimeInterval(msg)

def timeit(fn):
    def f(*args, **kwargs):
        with timing(f"{fn.__name__}"):
            return fn(*args, **kwargs)
    return f


# Global variable to hold the default temporary directory
DEFAULT_TEMP_DIR = None

def set_default_temp_dir(temp_dir):
    """
    Sets the default temporary directory.

    Args:
        temp_dir (str): Path to the default temporary directory.
    """
    global DEFAULT_TEMP_DIR
    DEFAULT_TEMP_DIR = temp_dir


def create_temp_file(write_data=None, suffix='', mode="w", temp_dir=None):
    """
    Creates a temporary file with optional data written to it.

    Args:
        write_data (str): Data to be written to the temporary file.
        suffix (str): Suffix for the temporary file.
        mode (str): File opening mode ('w' for text, 'wb' for binary).
        temp_dir (str): Directory to create the temporary file in. If None, uses the system default.

    Returns:
        _io.TextIOWrapper or _io.BufferedWriter: A file object for the temporary file.
    """
    # Create a temporary directory if not specified
    if temp_dir is None:
        temp_dir = DEFAULT_TEMP_DIR if DEFAULT_TEMP_DIR is not None else tempfile.gettempdir()


    # Ensure the temporary directory exists
    os.makedirs(temp_dir, exist_ok=True)

    # Create a named temporary file in the specified directory
    temp_file = tempfile.NamedTemporaryFile(suffix=suffix, delete=False, mode=mode, dir=temp_dir)

    if write_data:
        temp_file.write(write_data.encode() if mode == "w" else write_data)
        temp_file.flush()

    # It's important to return the file name for Nextflow to be able to access it
    return temp_file


def fasta_from_seq(name, seq):
    ## Input: Name, seq can be str or iterable yielding str
    result = []
    if not isinstance(name, str):
        try:
            for n, s in zip(name, seq):
                result.append('>{}\n{}'.format(n, s))
        except TypeError:
            logger.error('Please provide a iterable or string')
            raise TypeError
    else:
        result = ['>{}\n{}'.format(name, seq)]
    return '\n'.join(result)+'\n'
