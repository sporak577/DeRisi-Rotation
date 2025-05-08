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
    - percent recovered: length of the reassembled fragment relative to the full original protein
    - percent overlap recovered: propoertion of the original protein sequence covered by the fragment
    - percent recovered when summing over all fragments of one protein relative to the full original protein.

Output:
- a FASTA file of all reassembled fragments 
- a CSV summary file with one row per fragment, reporting identity and recovery statistics .
    - Fragment length 
    - Alignment identity
    - Recovery statistics
    - Cumulative coverage of the full protein
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re
from collections import defaultdict
import os

from Bio import Align
aligner = Align.PairwiseAligner()
aligner.mode = 'global'  # mimics globalxx behavior
aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = 0
aligner.extend_gap_score = 0

from itertools import chain

date = "050725"

# ----- INPUT FILES ------
original_proteins_fasta = '/Users/sophieporak/Library/CloudStorage/Box-Box/DeRisi/0.94_final_library_protein_records.faa' #Full-length original protein sequences
input_nt_tiles = "/Users/sophieporak/Library/CloudStorage/Box-Box/DeRisi/FINAL LIBRARY ORDERED/0.96_tiling_out_nt_tiles_cp_cleaned_no_linkers.fasta"
output_dir = f'/Users/sophieporak/Desktop/reassembly_output_{date}'

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
    header_main = record.id  # e.g. XRB79013.1_35_1
    parts = header_main.split("_")

    # Extract tile number safely
    try:
        tile_num = int(parts[-2])  # second-to-last part is always the tile number
    except ValueError:
        print(f"Could not extract tile number from: {record.id}")
        continue

    # Extract protein ID by joining everything before the tile number
    protein_id = "_".join(parts[:-2])
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
            SeqRecord(Seq(frag), id=frag_id, description=f"from {protein_id}")
        )

# save reassembled proteins
os.makedirs(os.path.dirname(reassembled_fasta), exist_ok=True)

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

    # Align, with mismatch and gap = 0 
    alignments = aligner.align(orig_seq, re_seq)
    alignment = alignments[0] #regions on the original sequence

    # Percent identity = exact matches / length of aligned region in fragment
    identity = alignment.score / len(re_seq) * 100

    # Extract alignment span on the original sequence
    orig_spans = alignment.aligned[0]  # List of (start, end) tuples on original. each list contains (start, end) positions of aligned blocks (gaps can split them)
    alignment_range = sum(end - start for start, end in orig_spans) #this gives me the total sum of the blocks that are aligned
    overlap_recovered = alignment_range / len(orig_seq) * 100

    # Recovery = length of fragment / length of full protein
    recovered = len(re_seq) / len(orig_seq) * 100
    

    rows.append({
        "FragmentID": frag_id,
        "ProteinID": protein_id,
        "OriginalLength": len(orig_seq),
        "FragmentLength": len(re_seq),
        "PercentIdentity": round(identity, 2), #number of matches / fragment length
        "PercentRecovered": round(recovered, 2), #how much did I reconstruct?
        "PercentOverlapRecovered": round(overlap_recovered, 2) #how much of the original protein did this fragment align to?
    })

# ------ STEP 4b: CALCULATE CUMULATIVE COVERAGE PER PROTEIN --------
coverage_map = defaultdict(list)

# collect aligned spans for each protein
for record in SeqIO.parse(reassembled_fasta, "fasta"):
    frag_id = record.id
    protein_id = frag_id.split("_frag")[0]

    if protein_id not in originals: 
        continue 

    re_seq = str(record.seq)
    orig_seq = originals[protein_id]
    alignments = aligner.align(orig_seq, re_seq)
    alignment = alignments[0]

    for start, end in alignment.aligned[0]:
        coverage_map[protein_id].append((start, end))

# merge these intervals and compute coverage 
def merge_intervals(intervals):
    sorted_intervals = sorted(intervals)
    merged = []
    for current in sorted_intervals:
        if not merged or current[0] > merged[-1][1]: #this checks if the current interval overlaps or touches the previous one, stored in merged[-1]. 
            #if the start of the current interval is after the end of the last one in merged, there's no overlap so new block. 
            merged.append(list(current))
        else: 
            merged[-1][1] = max(merged[-1][1], current[1])
    return merged 

cumulative_coverage = {}
for prot_id, intervals in coverage_map.items():
    merged = merge_intervals(intervals)
    total_covered = sum(end - start for start, end in merged)
    prot_len = len(originals[prot_id])
    cumulative_coverage[prot_id] = round(total_covered / prot_len * 100, 2)

for row in rows:
    row["CumulativeCoverage"] = cumulative_coverage.get(row["ProteinID"], 0.0)

# ----- STEP 5: EXPORT THE CSV SUMMARY -----
df = pd.DataFrame(rows)
df = df.sort_values(by="CumulativeCoverage", ascending=False) #sort by percent recovered, so essentially how much of the original protein got reconstructed
os.makedirs(output_dir, exist_ok=True)
df.to_csv(summary_output, index=False)
print(f"summary written to {summary_output}")



    
