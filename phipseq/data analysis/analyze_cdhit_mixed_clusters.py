"""
This script requires 

Input: 
- load FASTA file prior to clustering
- .clstr file 

Output: 
- file that checks for strain/family mixing upon clustering tiles
"""

import re 
from collections import defaultdict 
import csv
from Bio import Seq.IOError

# ----- Required Inputs -----
fasta_file = "your_input.fasta"
clstr_file = "your_file.clstr"
output_file = "clusters_with_mixed_strains.csv"

# ----- STEP 1: PARSE METADATA FROM FASTA HEADERS   
tile_to_metadata = {}

for record in SeqIO.parse(fasta_file, "fasta"):
    tile_id = record.id #e.g.WFG38034.1_1 from  >WFG38034.1_1 WFG38034.1 glycoprotein precursor|Aba-Mianyang virus|SC/C3-30.18/2021|Ochotona sp.|China|Aug-2021|S|?| | tile 1 of 23

    #Description: WFG38034.1 glycoprotein precursor|Aba-Mianyang virus|SC/C3-30.18/2021|Ochotona sp.|China|Aug-2021|S|?| | tile 1 of 23
    desc = record.description
    parts = desc.split("|")
    
    try: 
        id_and_name = parts[0].split()
        accession = id_and_name[0]
        protein_name = " ".join(id_and_name[1:]) #glycoprotein precursor 
        strain = parts[1].strip()
        isolate = parts[2].strip()
        host = parts[3].strip()
        country = parts[4].strip()
        collection_date = parts[5].strip()
    except IndexError:
        print(f"Skipping malformed header: {desc}")
        continue 

    tile_to_metadata[tile_id] = {
        "accession": accession, 
        "strain": strain, 
        "host": host,
        "country": country,
        "collection_date": collection_date
    }

# ------ STEP2: PARSE CD-HIT .clstr file ------
clusters = []
current_cluster = []

with open(clstr_file) as f:
    for line in f: 
        line = line.strip()
        if line.startswith(">Cluster"):
            if current_cluster: #checks if the list is non-empty, e.g. if current_cluster contains any elements do something. 
                clusters.append(current_cluster)
                current_cluster = []
        else:
            """
            (.*?) . matches any character (except newlines), * means "zero or more times", ? makes it non greedy - stops as soon as possible
            so this captures everything between > and ...
            \.\.\. matches a literal ... (each . must be escaped with \ in regex)
            """
            match = re.search(f">(.*?)\.\.\.", line) 
            if match: 
                tile_id = match.group(1) #0	46aa, >WFG38034.1_10... *, captures 'WFG38034.1_10'
                current_cluster.append(tile_id)
        
        if current_cluster:
            clusters.append(current_cluster)
        
# ------- STEP3: DETECT CLUSTERS WITH MIXED STRAINS/ACCESSIONS -------

cluster_analysis = []

for i, cluster in enumerate(clusters):
    strains = set()
    accessions = set() 
    hosts = set()
    countries = set()

    for tile in cluster: 
        meta = tile_to_metadata.get(tile)
        if not meta:
            continue
        strains.add(meta["strain"])
        accessions.add(meta["accession"])
        hosts.add(meta["host"])
        countries.add(meta["country"])

    #classify cluster type
    if len(accessions) == 1:
        category = "Same protein"
    elif len(strains) == 1:
        category = "Different proteins (same strain)"
    else: 
        category = "Different strains"

    cluster_analysis.append({
        "cluster_id": f"Cluster_{i}",
        "num_tiles": len(cluster),
        "cluster_type": category,
        "strains": "; ".join(sorted(strains)),
        "accessions": "; ".join(sorted(accessions)),
        "hosts": "; ".join(sorted(hosts)),
        "countries": "; ".join(sorted(countries)),
        "tiles": "; ".join(cluster)
    })

