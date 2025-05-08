
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

input_fasta = "/Users/sophieporak/Desktop/Malaria_MostSeroreactive_Haleigh.fa"     # Change to your input file path
output_fasta = "Malaria_MostSeroreactive_trimmed_46aa.fasta"

trimmed_records = []

for record in SeqIO.parse(input_fasta, "fasta"):
    seq = str(record.seq)
    if len(seq) < 46:
        print(f"Skipping {record.id}: sequence too short ({len(seq)} aa)")
        continue

    # Calculate how many residues to trim from each side
    excess = len(seq) - 46
    left_trim = excess // 2
    right_trim = excess - left_trim

    trimmed_seq = seq[left_trim:len(seq) - right_trim]

    new_record = SeqRecord(
        Seq(trimmed_seq),
        id=record.id,
        description=f"{record.description} | trimmed to 46aa"
    )
    trimmed_records.append(new_record)

# Write output
SeqIO.write(trimmed_records, output_fasta, "fasta")
print(f"Trimmed sequences saved to: {output_fasta}")
