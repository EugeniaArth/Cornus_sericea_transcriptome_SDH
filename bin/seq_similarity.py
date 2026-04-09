from Bio import SeqIO
from Bio import pairwise2

# check sequence similarity of DQD_SDH (nucleotide) 

# Load sequences
records = list(SeqIO.parse("/Files/DQD_SDH_possible_nt.fasta", "fasta"))
seqs = [str(r.seq) for r in records]
names = [r.id for r in records]

def percent_identity(seq1, seq2):
    aln = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]
    matches = sum(a == b for a, b in zip(aln.seqA, aln.seqB))
    return matches / len(aln.seqA) * 100

# Pairwise comparisons
for i in range(len(seqs)):
    for j in range(i+1, len(seqs)):
        pid = percent_identity(seqs[i], seqs[j])
        print(f"{names[i]} vs {names[j]}: {pid:.2f}%")
