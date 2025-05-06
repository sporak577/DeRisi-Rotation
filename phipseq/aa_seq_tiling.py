'''
This script will read in a fasta file of protein sequences and
1) deduplicate -> save in fasta file
2) tile them to 46aa, with 23ß aa overlap
2.2) deduplicate again on tile level & get rid of tiles containing X & other invalid characters -> save header and sequence in fasta file
3) codon optimize 
4) purge restriction sites by using synonymous mutations
5) write to output FASTA

it will keep metadata in FASTA headers of generated tiles

it will also write an output FASTA of generated peptides before peptide processing!
'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import pandas as pd
import numpy as np

import os 
import random

cd_hit_threshold = 0.96

fasta_file = '/Users/sophieporak/Library/CloudStorage/Box-Box/DeRisi/cd-hit_0.96_on_0.94_aa_tiles_final_library/cd-hit_0.96_0.94_aa_tiles.fasta' #update path 

# path to output directory
#output_dir = f"/Users/sophieporak/Documents/DeRisi_data /cd-hit arenavirus merged/{cd_hit_threshold}_cd-hit_tiles"


#fasta_file = "/Users/sophieporak/Documents/DeRisi_data /arenavirus_merged.fasta" #update path 

# path to output directory
output_dir = '/Users/sophieporak/Desktop/cd-hit_0.96_0.94_final_library_protein_records_tiling_out_050525'

# create directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# name of output files
output_tiling = os.path.join(output_dir, f'{cd_hit_threshold}_tiling_out_aa_preprocess_tiles.fasta')
output_fasta = os.path.join(output_dir, f'{cd_hit_threshold}_tiling_out_nt_tiles.fasta')
output_duplicate_proteins = os.path.join(output_dir,f'{cd_hit_threshold}_tiling_out_duplicate_proteins.fasta')
output_X_tiles = os.path.join(output_dir, f'{cd_hit_threshold}_tiles_aa_with_X_characters.fasta')
output_duplicate_tiles = os.path.join(output_dir, f'{cd_hit_threshold}_tiling_out_duplicate_aa_tiles.fasta')
translation_mismatched_tiles = os.path.join(output_dir, f'{cd_hit_threshold}_tiling_out_translation_mismatch.fasta')
output_unresolved_restriction = os.path.join(output_dir, f'{cd_hit_threshold}_unresolved_restriction_sites.fasta')



# ===== Codon Optimization to E.Coli, copied from Haleigh Miller ======

AA2NA = {
"A": list("GCT,GCC,GCA,GCG".split(",")),
"R": list("CGT,CGC,CGA,CGG,AGA,AGG".split(",")),
"N": list("AAT,AAC".split(",")),
"D": list("GAT,GAC".split(",")),
"C": list("TGT,TGC".split(",")),
"Q": list("CAA,CAG".split(",")),
"E": list("GAA,GAG".split(",")),
"G": list("GGT,GGC,GGA,GGG".split(",")),
"H": list("CAT,CAC".split(",")),
"I": list("ATT,ATC,ATA".split(",")),
"L": list("TTA,TTG,CTT,CTC,CTA,CTG".split(",")),
"K": list("AAA,AAG".split(",")),
"M": list("ATG".split(",")),
"F": list("TTT,TTC".split(",")),
"P": list("CCT,CCC,CCA,CCG".split(",")),
"S": list("TCT,TCC,TCA,TCG,AGT,AGC".split(",")),
"T": list("ACT,ACC,ACA,ACG".split(",")),
"W": list("TGG".split(",")),
"Y": list("TAT,TAC".split(",")),
"V": list("GTT,GTC,GTA,GTG".split(",")),
"*": list("TAA,TGA,TAG".split(","))}


NA2AA = {'GCT': 'A',
 'GCC': 'A',
 'GCA': 'A',
 'GCG': 'A',
 'CGT': 'R',
 'CGC': 'R',
 'CGA': 'R',
 'CGG': 'R',
 'AGA': 'R',
 'AGG': 'R',
 'AAT': 'N',
 'AAC': 'N',
 'GAT': 'D',
 'GAC': 'D',
 'TGT': 'C',
 'TGC': 'C',
 'CAA': 'Q',
 'CAG': 'Q',
 'GAA': 'E',
 'GAG': 'E',
 'GGT': 'G',
 'GGC': 'G',
 'GGA': 'G',
 'GGG': 'G',
 'CAT': 'H',
 'CAC': 'H',
 'ATT': 'I',
 'ATC': 'I',
 'ATA': 'I',
 'TTA': 'L',
 'TTG': 'L',
 'CTT': 'L',
 'CTC': 'L',
 'CTA': 'L',
 'CTG': 'L',
 'AAA': 'K',
 'AAG': 'K',
 'ATG': 'M',
 'TTT': 'F',
 'TTC': 'F',
 'CCT': 'P',
 'CCC': 'P',
 'CCA': 'P',
 'CCG': 'P',
 'TCT': 'S',
 'TCC': 'S',
 'TCA': 'S',
 'TCG': 'S',
 'AGT': 'S',
 'AGC': 'S',
 'ACT': 'T',
 'ACC': 'T',
 'ACA': 'T',
 'ACG': 'T',
 'TGG': 'W',
 'TAT': 'Y',
 'TAC': 'Y',
 'GTT': 'V',
 'GTC': 'V',
 'GTA': 'V',
 'GTG': 'V',
 'TAA': '*',
 'TGA': '*',
 'TAG': '*'}

def tiling(seq, tile_len=46, overlap=23): 
    """
    Tiles the input amino acid sequence into overlapping peptides. 
    
    Parameters
    ----------
    seq : str
        amino acid sequence 
    tile_len : int
        length of each tile (default 48)
    overlap : int
        number of overlapping amino acids between tiles (default 24)
    
    Returns
    -------
    tiles : list of str
        list of overlapping amino acid tiles 
    """
    step = tile_len - overlap
    tiles = []
    for i in range(0, len(seq) - tile_len + 1, step):
        tiles.append(seq[i:i + tile_len])

    #add the final tile to include any trailing residues 
    if len(seq) > 0 and (len(seq) - tile_len) % step != 0: 
        final_tile = seq[-tile_len:]
        if final_tile not in tiles: 
            tiles.append(final_tile)
    
    return tiles




def aa2na(seq):
    """
    kept from Haleigh Miller 

    translates amino acid to e.coli preferred codons, randomly choosing codon for each AA
    
    Parameters
    ----------
    seq: str
        amino acid sequence
    
    Returns
    -------
    seq: str
        nucleotide sequence
    """
    na_seq = [random.choice(AA2NA.get(c, ["---"])) for c in seq]
    return "".join(na_seq)

            
def replace_restriction_sites(seq):
    """
    Iteratively replaces all restriction enzyme cut sites using synonymous mutations (in-frame codon)
    
    Parameters
    ----------
    seq: str
        codon-optimized nucleotide sequence 

    AA2NA : dict
        amino acid to synonymous codon mapping
    
    NA2AA : dict
        codon to amino acid mapping
    
    Returns
    -------
    new_seq: str or None
        nucleotide sequence with synonomous mutation at restriction site, or None if no restriction sites
    
    Notes
    -----
    restriction sites:
    ecoRI='GAATTC'
    hindIII='AAGCTT'
    bamHI='GGATCC'
    XhoI='CTCGAG'
    """
    restriction_sites=['GAATTC','AAGCTT','GGATCC','CTCGAG']
    max_passes = 100 
    
    for _ in range(max_passes):
        #checks if we successfully removed any site in this round, if not stop early.
        modified = False

        for site in restriction_sites: 
            #finds the first index where the site appears
            idx = seq.find(site)
            #if not found, skip and go to next one
            if idx == -1: 
                continue 

            #try mutating one overlapping codon within the site window
            for offset in range(len(site) -2): #slide window to find codon overlap
                #this snaps to the nearest codon boundary (multiple of 3)
                codon_start = (idx + offset) - ((idx + offset) % 3)
                #exract a 3-letter codon starting at that position
                codon = seq[codon_start:codon_start + 3]

                #validate codon
                if len(codon) != 3 or codon not in NA2AA: 
                    continue #malformed codon or not standard, invalid or isn't in reverse mapping so skip it. 

                #find the amino acid that the codon encodes
                aa = NA2AA[codon]
                #get all synonymous codons for that amino acid
                synonyms = AA2NA[aa]

                #skip if there are no alternatives e.g. TGG -> Tryptophan
                if len(synonyms) <= 1:
                    continue #no alternative codons for this AA

                #try each alternative codon and construct a sequence with the substituted codon
                for alt in synonyms: 
                    if alt != codon: 
                        #mutate codon and check if site is broken
                        mutated_seq = seq[:codon_start] + alt + seq[codon_start + 3:]
                        if site not in mutated_seq:
                            seq = mutated_seq
                            modified = True 
                            break 

                    if modified: 
                        break #codon successfully mutated 
                if modified: 
                    break #move on to re-scan all sites 
            if not modified:
                break

    return seq
            

# ===== Main Loop =====

seen_proteins = set()
records_out = []
aa_records_out = []
seen_peptides = set()
duplicate_count = 0

# full protein duplicates out
full_protein_duplicates_out = []

# detect tiles with X and invalid characters 
tiles_with_x_out = []
# detect duplicate tiles 
duplicate_tiles_out = []
#translation mismatches 
translation_mismatches = []
unresolved_restriction_tiles = []


for record in SeqIO.parse(fasta_file, "fasta"):
    aa_seq = str(record.seq)

    # skip duplicates early
    if aa_seq in seen_proteins: 
        full_protein_duplicates_out.append(record) #save for output
        continue  #skips the rest of this loop iteration (don't tile or process) and move to the next one
    seen_proteins.add(aa_seq)

    tiled_peptides = tiling(aa_seq)


    for i, peptide in enumerate(tiled_peptides):
        # Fasta header ID
        header = f"{record.id}_{i+1}"
        # Metadata in the record.description
        desc = f"{record.description} | tile {i+1} of {len(tiled_peptides)}"

        valid_aa = set("ARNDCQEGHILKMFPSTWYV")
        if not set(peptide).issubset(valid_aa) or "X" in peptide: #is every element in set peptide also in set valid_aa?
            tiles_with_x_out.append(SeqRecord(Seq(peptide), id=header, description = desc))
            continue

        if peptide in seen_peptides: 
            duplicate_count += 1
            duplicate_tiles_out.append(SeqRecord(Seq(peptide), id=header, description = desc))
            continue 

        seen_peptides.add(peptide)

        na_seq = aa2na(peptide) #codon optimization
        na_seq_clean = replace_restriction_sites(na_seq) or na_seq

        # Final check for any residual restriction site
        restriction_sites = ['GAATTC', 'AAGCTT', 'GGATCC', 'CTCGAG']
        if any(site in na_seq_clean for site in restriction_sites):
            print(f"Unresolved restriction site in: {header}")
            unresolved_restriction_tiles.append(SeqRecord(Seq(peptide), id=header, description=desc + " | unresolved restriction site"))
            continue  # skip this tile!

        #sanity check: translated sequences match original peptide
        translated = str(Seq(na_seq_clean).translate())
        if translated != peptide:
            translation_mismatches.append((header, peptide, translated))
            print(f"Translation mismatch for {header} - mutation may have broken reading frame or cause mis-translation")
            print(f"original: {peptide}")
            print(f"translated: {translated}")
            #skip saving this tile!
            continue 

        # Create SeqRecord, wraps the nucleotide sequence string into a Seq object. id = header means becomes the identifier in the FASTA, the part right after >. 
        #description=desc becomes the rest of the FASTA header line, holding metadata like protein name, virus, location etc. 
        rec = SeqRecord(Seq(na_seq_clean), id=header, description=desc)
        # Create amino acid FASTA entry (for validation), unprocessed
        aa_rec = SeqRecord(Seq(peptide), id=header, description=desc) 

        records_out.append(rec)
        aa_records_out.append(aa_rec)

        

#backtranslation sanity check
output_backtranslated = os.path.join(output_dir, f'{cd_hit_threshold}_tiling_out_backtranslated_aa.fasta')
backtranslated_records = []

for rec in records_out:
    back_aa = str(rec.seq.translate())
    aa_rec = SeqRecord(Seq(back_aa), id=rec.id, description=rec.description + " | backtranslated")
    backtranslated_records.append(aa_rec)

with open(output_backtranslated, "w") as out_back:
    SeqIO.write(backtranslated_records, out_back, "fasta")

print(f"Done! {len(backtranslated_records)} backtranslated AA tiles written to '{output_backtranslated}'")

with open(output_fasta, "w") as out_f:
    SeqIO.write(records_out, out_f, "fasta")

print(f"Done! {len(records_out)} unique tiles written to '{output_fasta}'")

with open(output_tiling, "w") as out_aa:
    SeqIO.write(aa_records_out, out_aa, "fasta")

print(f"Done! {len(aa_records_out)} unique tiles written to '{output_tiling}'")

with open(output_duplicate_proteins, "w") as out_dup_prots:
    SeqIO.write(full_protein_duplicates_out, out_dup_prots, "fasta")

print(f"Skipped {len(full_protein_duplicates_out)} duplicate full-length proteins written to '{output_duplicate_proteins}'")

with open(output_X_tiles, "w") as out_x:
    SeqIO.write(tiles_with_x_out, out_x, "fasta")

print(f"Excluded {len(tiles_with_x_out)} tiles with X character in them. Written to '{output_X_tiles}'")

with open(output_duplicate_tiles, "w") as out_dup_tiles:
    SeqIO.write(duplicate_tiles_out, out_dup_tiles, "fasta")

print(f"Skipped {len(duplicate_tiles_out)} duplicate tiles written to '{output_duplicate_tiles}'")


with open(translation_mismatched_tiles, "w") as out_mismatch_tiles: 
    SeqIO.write(translation_mismatches, out_mismatch_tiles, "fasta")
print(f"Skipped {len(translation_mismatches)} translation mismatched tiles written to '{output_duplicate_tiles}'")

with open(output_unresolved_restriction, "w") as out_unresolved:
    SeqIO.write(unresolved_restriction_tiles, out_unresolved, "fasta")
print(f"Excluded {len(unresolved_restriction_tiles)} tiles due to unresolved restriction sites. Written to '{output_unresolved_restriction}'")

