"""
Require inputs:
- original proteins FASTA (pre-tiling)
- FASTA of nucleotide-tiles (linkers removed)
- need to have shared protein IDs in headers 

Process:
1) extract protein_id, tile_number and nt sequence
2) translate nt sequence to peptide
3) sort by tile number, then stitch them together based on overlap that you previously defined
4) percent identity, so the amino acid match over aligned region
5) percent recovery (reassembled length / original length)
6) list of mismatches or missing regions

Output:
- a csv summary table
"""

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re
from collections import defaultdict
import os

# ----- INPUT FILES ------
original_proteins_fasta = 'original_proteins.fasta' #Full-length original protein sequences
input_nt_tiles = "nt_tiles_no_linkers.fasta"
output_dir = "reassembly_output"

# ----- OUTPUT FILES -----
reassembled_fasta = os.path.join(output_dir, "reassembled_proteins.fasta")
summary_output = os.path.join(output_dir, "reassembly_comparison_summary.csv")

# ----- PARAMETERS -------

tile_len = 46 
overlap = 23 


# ----- STEP 1: PARSE AND TRANSLATE TILES, GROUP BY PROTEIN ID ------
protein_tiles = defaultdict(list)

for record in SeqIO.parse(input_nt_tiles, "fasta"):
    desc = record.description
    # exctract second column from e.g. >XRB79013.1_35_1 XRB79013.1_35 XRB79013.1 RNA-dependent RNA polymerase|Mammarenavirus wenzhouense|H2024-1|Mus musculus|China|2024|?|?| | tile 35 of 96 | tile 1 of 1_trimmed
    fields = record.description.split()
    protein_id = fields[0].split("_")[0]
    tile_num = int(fields[1].split("_")[1])  # Extract tile number from 2nd field
    peptide = str(Seq(record.seq).translate())

    protein_tiles[protein_id].append((int(tile_num), peptide))

# ----- STEP 2: REASSEMBLE PEPTIDE SEQUENCES FROM TILES -----
reassembled_records = []

for protein_id, tiles in protein_tiles.items():
    sorted_tiles = sorted(tiles, key=lambda x:x[0])

    fragments = []
    current_seq = sorted_tiles[0][1]

    for i in range(1, len(sorted_tiles)):
        prev = current_seq[-overlap:]
        curr = sorted_tiles[i][1]
        if prev == curr[:overlap]:
            current_seq += curr[overlap:] # continue building fragment
        else:
            fragments.append(current_seq) # save current (valid so far) fragment
            print(f"Fragment break in {protein_id} at tile {i+1}")
            current_seq = curr  # start new fragment
    
    fragments.append(current_seq)  # add the last one



    
