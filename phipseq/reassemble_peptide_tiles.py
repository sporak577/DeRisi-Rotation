"""
This script evaluates how well tiled nucleotide sequences (without linkers) reconstruct

Required inputs:
- original proteins FASTA (pre-tiling)
- FASTA of nucleotide-tiles (linkers removed)
- need to have shared protein IDs in headers 

Process:
1) extract protein_id, tile_number and nt sequence
2) translate nucleotide sequence into peptide sequences
3) sort by tile number, then stitch them together based on overlap that you previously defined
  - if two adjacent tiles do not share a perfect overlap, a new fragment is started 
4) for each reassembled fragment:
  - align it to the corresponding full-length protein.
  - calculate the following metrics: 
    - percent identity: proportion of exact amino acid matches in the alignment. 
    - pecent recovered: length of the reassembled fragment relative to the full original protein
    - percent overlap recovered: propoertion of the original protein sequence covered by the fragment

Output:
- a FASTA file of all reassembled fragments 
- a CSV summary file with one row per fragment, reporting identity and recovery statistics .
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

    #save each fragment separately with index
    for j, frag in enumerate(fragments):
        frag_id = f"{protein_id}_frag{j+1}"
        reassembled_records.append(
            SeqRecord(Seq(frag), id=frag_id, description=f"Fragment {j+1} from {protein_id}")
        )

# save reassembled proteins
with open(reassembled_fasta, "w") as out_f:
    SeqIO.write(reassembled_records, out_f, "fasta")

# ----- STEP 3: LOAD ORIGINAL PROTEINS FOR COMPARISON -----
originals = {record.id.split()[0]: str(record.seq) for record in SeqIO.parse(original_proteins_fasta, "fasta")}

# ----- STEP4: COMPARE ORIGINAL TO REASSEMBLED -----
rows = []
for record in SeqIO.parse(reassembled_fasta, "fasta"):
    frag_id = record.id 
    protein_id = frag_id.split("_frag")[0] # base ID

    if protein_id not in originals:
        print(f"original not found for fragment {frag_id}")
        continue

    re_seq = str(record.seq)
    orig_seq = originals[protein_id]

    # Align with globalxx (identity-based)
    alignment = pairwise2.align.globalxx(orig_seq, re_seq, one_alignment_only=True)[0]
    alignment_range = alignment[4] - alignment[3] #end - begin
    overlap_recovered = alignment_range / len(orig_seq) * 100 
    identity = alignment[2] / alignment[4] * 100 #matches / aligned length, number of exact amino acid matches over the length of alignment
    recovered = len(re_seq) / len(orig_seq) * 100 #length of reassembled fragment vs original

    rows.append({
        "FragmentID": frag_id,
        "ProteinID": protein_id,
        "OriginalLength": len(orig_seq),
        "FragmentLength": len(re_seq),
        "PercentIdentity": round(identity, 2), 
        "PercentRecovered": round(recovered, 2), #how much did I reconstruct?
        "PercentOverlapRecovered": round(overlap_recovered, 2) #how much of the original protein did this fragment align to?
    })

# ----- STEP 5: EXPORT THE CSV SUMMARY -----
df = pd.DataFrame(rows)
df = df.sort_values(by="PercentRecovered", ascending=False) #sort by percent recovered, so essentially how much of the original protein got reconstructed
os.makedirs(output_dir, exist_ok=True)
df.to_csv(summary_output, index=False)
print(f"summary written to {summary_output}")



    
