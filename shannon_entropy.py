"""
Authors: Sophie-Christine Porak & ChatGPT
This script will assess how peptides score in the following domains: 
1) Diversity - Shannon Entropy

Input: 
FASTA file 

Output: 
 - CSV file: Shannon entropy distribution of sequences ranked lowest to highest entropy.
 - Histogram: Shannon entropy distribution of sequences. 

"""
from Bio import SeqIO 
from collections import Counter
import math 
import csv
import matplotlib.pyplot as plt

input_file = "test_entropy.fasta"
looking_at_which_viruses = "test"
looking_at_nt_or_aa = "nucleotides"
date = "05-02-25"

def shannon_entropy(seq):
    #Counter(seq) counts how many times each amino acid or nucleotide appears. 
    counts = Counter(seq)
    #get total length of sequence
    total = len(seq)
    """
    for each amino acid it takes the frequency (c / total) and then multiplies by log2(frequency). this is being summed over all amino acids.
    some background: in information theory the entropy of a random variable quantifies the average level of uncertainty with the variable's potential states.
    so the more diverse our sequence is, the higher the entropy. 
    calculated with: 

    H(X) = - sum[possibility of observing amino acid x * log2 (possibility of observing amino acid x)]
      where x is element of X, and X = all 20 amino acids.
    
    H(X) is in bits. For 20 equally probable amino acids Hmax is around 4.32 bits, which would be the maximal entropy in this scenario. 

    for nucleotides the maximum entropy is 2.0 bits, because there are 4 bases. 
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
plt.hist([x[1] for x in entropy_values], bins=30, edgecolor='black', alpha=0.8)


plt.title(f"Shannon Entropy Distribution of {looking_at_nt_or_aa} tiles, {looking_at_which_viruses}")
plt.xlabel("Shannon Entropy (bits)")
plt.ylabel("Number of Sequences")
plt.grid(True)
plt.tight_layout()
plt.savefig(f"entropy_distribution_{looking_at_which_viruses}.png", dpi=300)
plt.show()

#sort results 
entropy_values.sort(key=lambda x: x[1], reverse=False)

#write entropy-ranked output
with open(f"entropy_rank_{looking_at_which_viruses}_{looking_at_nt_or_aa}_{date}.csv", "w", newline="") as out_csv:
    writer = csv.writer(out_csv)
    writer.writerow(["Record ID", "Shannon Entropy", "Sequence"])
    writer.writerows(entropy_values)

print("Output files written")