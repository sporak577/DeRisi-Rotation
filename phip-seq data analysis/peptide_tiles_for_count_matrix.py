"""
Authors: ChatGPT & Sophie-Christine Porak 
"""
import csv
from Bio.Seq import Seq

input_fasta = '/Users/sophieporak/Library/CloudStorage/Box-Box/DeRisi/lassa fever project/FINAL LIBRARY ORDERED/0.96_tiling_out_nt_tiles_cp_cleaned_no_linkers.fasta'
output_csv = "parsed_tiles_with_translation.csv"

records = []
index = 0

with open(input_fasta) as f:
    seq_id = ref_seq_id = sequence = protein_description = virus_name = None

    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seq_id and sequence:
                # Translate nucleotide sequence to protein
                aa_seq = str(Seq(sequence).translate(to_stop=True))
                records.append([index, seq_id, ref_seq_id, protein_description, virus_name, aa_seq])
                index += 1
            parts = line[1:].split()
            seq_id = parts[0]
            ref_seq_id = parts[2]

            # Rebuild the description and parse fields
            description = " ".join(parts[3:])
            fields = description.split("|")
            protein_description = fields[0].strip() if len(fields) > 0 else "" #removing any leading and trailing whitespace
            virus_name = fields[1].strip() if len(fields) > 1 else ""
            sequence = ""
        else:
            sequence += line

    # Add the last record
    if seq_id and sequence:
        aa_seq = str(Seq(sequence).translate(to_stop=True))
        records.append([index, seq_id, ref_seq_id, protein_description, virus_name, aa_seq])

# Write to CSV
with open(output_csv, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["index", "seq_id", "ref_seq_id", "protein_description", "virus_name", "amino_acid_sequence"])
    writer.writerows(records)

print(f"Saved to {output_csv}")
