"""
Authors: ChatGPT & Sophie-Christine Porak 
"""
import csv
from Bio.Seq import Seq

input_fasta = 'your_file'
output_csv = "parsed_tiles_with_translation.csv"

records = []
index = 0

with open(input_fasta) as f:
    seq_id = ref_seq_id = sequence = protein_description = virus_name = None

    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seq_id and sequence:
                aa_seq = str(Seq(sequence).translate(to_stop=True))
                fragment_number = seq_id.split("_")[1] if len(seq_id.split("_")) > 1 else ""
                records.append([index, seq_id, ref_seq_id, fragment_number, protein_description, virus_name, aa_seq])
                index += 1

            parts = line[1:].split()
            seq_id = parts[0]         # e.g. WFG38034.1_3_1
            ref_seq_id = parts[2]     # e.g. WFG38034.1

            # Extract description and virus info
            description = " ".join(parts[3:])
            fields = description.split("|")
            protein_description = fields[0].strip() if len(fields) > 0 else ""
            virus_name = fields[1].strip() if len(fields) > 1 else ""
            sequence = ""
        else:
            sequence += line

    # Final record
    if seq_id and sequence:
        aa_seq = str(Seq(sequence).translate(to_stop=True))
        fragment_number = seq_id.split("_")[1] if len(seq_id.split("_")) > 1 else ""
        records.append([index, seq_id, ref_seq_id, fragment_number, protein_description, virus_name, aa_seq])

# Write CSV
with open(output_csv, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["index", "seq_id", "ref_seq_id", "fragment_number", "protein_description", "virus_name", "amino_acid_sequence"])
    writer.writerows(records)

print(f"Saved to {output_csv}")
