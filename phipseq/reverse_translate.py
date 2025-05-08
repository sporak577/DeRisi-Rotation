"""
this is just a dummy check to translate your generated nucleotide tiles 
"""

from Bio import SeqIO
from Bio.Seq import Seq 

#here checking for WFG38034.1, and WFG38034.1

#paste nucleotide sequences here
nt_sqs = ["GTGGTTGGTGCTGTAGGAGCAATGTCCAGGGTATGCCTATGTCTCCTCTTCTCCGGCCTCTTATTATGGCGCGCGGCCGAACTGAGGAACTTGATTGAGCTGAAAATTGAGTGCCCGCATACCATTGGGCTAGGACAAGGATTGGTTATTGGATCTGTGTGATAAGCATATGCCATGGCCTC"]
for i, nt in enumerate(nt_sqs, start=1):
    #this stops translating at the first stop codon
    aa_seq = Seq(nt).translate(to_stop=True)
    print(f"Sequence {i}:")
    print(f" Nucleotide: {nt}")
    print(f" Amino acid: {aa_seq}\n")