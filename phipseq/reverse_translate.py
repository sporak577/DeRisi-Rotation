"""
this is just a dummy check to translate your generated nucleotide tiles 
"""

from Bio import SeqIO
from Bio.Record import SeqRecord
from Bio.Seq import Seq 


#paste nucleotide sequences here
nt_sqs = [ "",
         ""

]

for i, nt in enumerate(nt_sqs, start=1):
    #this stops translating at the first stop codon
    aa_seq = Seq(nt).translate(to_stop=True)
    print(f"Sequence {i}:")
    print(f" Nucleotide: {nt}")
    print(f" Amino acid: {aa_seq}\n")