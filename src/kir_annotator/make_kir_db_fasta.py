from kir_annotator.common import fasta_from_seq
import dill as pickle

db_path = 'data/kir.pickle'
output_path = 'data/kir_db.fasta'

with open(db_path, 'rb') as f:
    anno = pickle.load(f)
    
seqs = []
for g in anno[0].values():
    for a in g.alleles.values():
        seqs.append((g.name+'*'+a.name, a.seq))
        
with open(output_path, 'w') as f:
    f.write(fasta_from_seq(*zip(*seqs)))