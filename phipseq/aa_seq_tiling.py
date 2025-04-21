'''
This script will read in a fasta file of protein sequences and
1) tile them to 48aa, with 24 aa overlap
2) collapses on sequence identity
3) codon optimize 
4) purge restriction sites by using synonymous mutations
'''

from Bio import SeqIO
import pandas as pd
import numpy as np
import random


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

def tiling(seq, tile_len=48, overlap=24): 
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








def collapse():



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
    kept from Haleigh Miller 

    checks sequence for restriction enzyme cut sites and replaces first in-frame codon with synonomous mutation.
    
    Parameters
    ----------
    seq: str
        nucleotide sequence 
    
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
    
    new_seq=None
    for r in restriction_sites: 
        if seq.find(r) != -1:
            ##find codon at restriction site
            x=seq.find(r)
            n=x%3 #find start of codon
            codon=seq[x-n:x-n+3]
            
            ##replace with synonomous mutation
            tmp=AA2NA[NA2AA[codon]][:] #get copy of codon list
            if len(tmp)>1:
                tmp.remove(codon) #remove the one causing restriction site
                new_codon=tmp[0] #replace 
                new_seq=seq[:x-n]+new_codon+seq[x-n+3:]
                
            #tryptophan has only one codon so you cant do this
            else: 
                codon=seq[x-n+3:x-n+6]
                print(codon)
                tmp=AA2NA[NA2AA[codon]][:]
                tmp.remove(codon) #remove the one causing restriction site
                new_codon=tmp[0] #replace 
                new_seq=seq[:x-n]+new_codon+seq[x-n+3:]
                
            break
        
        else:
            continue 
    
    if new_seq:
        return new_seq
    else:
        return None
            

   