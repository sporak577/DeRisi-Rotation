"""
Authors: Sophie-Christine Porak & ChatGPT
This script will assess how peptides score in the following domains: 
1) Diversity - Shannon Entropy

Input: 
FASTA file 

Output: 

"""
from Bio import SeqIO 
from collections import Counter
import math 
import csv
import matplotlib.pyplot as plt

input_file = "/Users/sophieporak/Library/CloudStorage/Box-Box/DeRisi/Arenavirus/20250428/0.94_txid11617[Organism:exp]_042528_tiling_out/0.94_tiling_out_nt_tiles.fasta"
looking_at = "arenaviridae"

def shannon_entropy(seq):
    #Counter(seq) counts how many times each amino acid appear
    counts = Counter(seq)
    #get total length of sequence, here peptide
    total = len(seq)
    """
    for each amino acid it takes the frequency (c / total) and then multiplies by log2(frequency). this is being summed over all amino acids.
    some background: in information theory the entropy of a random variable quantifies the average level of uncertainty with the variable's potential states.
    so the more diverse our sequence is, the higher the entropy. 
    calculated with: 

    H(X) = - sum[possibility of observing amino acid x * log2 (possibility of observing amino acid x)]
      where x is element of X, and X = all 20 amino acids.
    
    H(X) is between 0 and 1. 

    note that the entropy is maximal if possibility of observing amino axid x is 1/20th for all amino acids, as this sums up to 1!
    """

    entropy = -sum((c/total) * math.log2(c/total) for c in counts.values())
    return entropy 

#store all entropy values
entropy_values = []

for record in SeqIO.parse(input_file, "fasta"):
    seq = str(record.seq)
    if len(seq) > 0:
        entropy = shannon_entropy(seq)
        entropy_values.append((record.id, round(entropy, 4), seq))

# Plot histogram 
plt.figure(figsize=(10, 6))
plt.hist(entropy_values, bins=30, edgecolor='black', alpha=0.8)
plt.title(f"Shannon Entropy Distribution of Peptide Sequences, {looking_at}")
plt.xlabel("Shannon Entropy")
plt.ylabel("Number of Sequences")
plt.grid(True)
plt.tight_layout()
plt.savefig(f"entropy_distribution_{looking_at}.png", dpi=300)
plt.show()

#sort results 
entropy_values.sort(key=lambda x: x[1], reverse=True)

#write entropy-ranked output
with open("entropy_ranked.csv", "w", newline="") as out_csv:
    writer = csv.writer(out_csv)
    writer.writerow(["Record ID", "Shannon Entropy", "Sequence"])
    writer.writerows(entropy_values)

print("Output files written")