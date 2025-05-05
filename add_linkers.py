"""
This script adds linkers to generated tiles (in nucleotide form)
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Define flanking sequences (as DNA)
flank_5 = "GTGGTTGGTGCTGTAGGAGCA"
flank_3 = "TGATAAGCATATGCCATGGCCTC"

# Input and output files
input_fasta = "/Users/sophieporak/Desktop/0.94_final_library_protein_records_tiling_out/0.94_tiling_out_nt_tiles.fasta"
output_fasta = "0.94_tiling_out_nt_tiles_with_flanks.fasta"

# Read input FASTA and add flanks
updated_records = []
for record in SeqIO.parse(input_fasta, "fasta"):
    new_seq = Seq(flank_5) + record.seq + Seq(flank_3)
    new_record = SeqRecord(new_seq, id=record.id, description="with 5' and 3' flanks")
    updated_records.append(new_record)

# Write to output FASTA
with open(output_fasta, "w") as out_handle:
    SeqIO.write(updated_records, out_handle, "fasta")

