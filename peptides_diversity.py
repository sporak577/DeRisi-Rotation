"""
Authors: Sophie-Christine Porak & ChatGPT
This script will assess how peptides score in the following domains: 
1) Repeated amino acid stretch length
2) Diversity - Shannon Entropy

Input: 
FASTA file 

Output: 

"""
from Bio import SeqIO 
from collections import Counter
import math 
import csv

input_file = "my_input.fasta"

def shannon_entropy(seq):
    counts = Counter(seq)
    total = len(seq)
    entropy = -sum((c/total) * math.log2(c/total) for c in counts.values())
    return entropy 

diversity_scores = []
repeat_scores = []

for record in SeqIO.parse(input_file, "fasta"):
    seq = str(record.seq)

    #Shannon entropy 
    entropy = shannon_entropy(seq)
    diversity_scores.append((record.id, round(entropy, 4), seq))

    #longest homopolymer stretch 
    max_run = 1 
    run = 1 
    for i in range(1, len(seq)):
        if seq[i] == seq[i -1]:
            run += 1
            max_run = max(max_run, run)
        else: 
            run = 1
    repeat_scores.append((record.id, max_run, seq))

#sort results 
diversity_scores.sort(key=lambda x: x[1], reverse=True)
repeat_scores.sort(key=lambda x: x[1], reverse=True)

#write entropy-ranked output
with open("entropy_ranked.csv", "w", newline="") as out_csv:
    writer = csv.writer(out_csv)
    writer.writerow(["Record ID", "Shannon Entropy", "Sequence"])
    writer.writerows(diversity_scores)

#write homopolymer-ranked output 
with open("longest_repeat_to_shortest.csv", "w", newline="") as out_csv:
    writer = csv.writer(out_csv)
    writer.writerow(["Record ID", "Longest AA Repeat", "Sequence"])
    writer.writerows(repeat_scores)

print("Output files written")