#!/bin/env Python
import sys
from Bio import SeqIO
import pandas as pd

seqs: str = ""
look_in = sys.argv[1]

busco_dir: str = f"{look_in}/run_sordariomycetes_odb10"
table = pd.read_csv(f"{busco_dir}/full_table.tsv", sep="\t",  skiprows=2)
complete = table[table["Status"] == "Complete"]
for gene in complete.iterrows():
    id = gene[1]['# Busco id']
    strand = gene[1]["Strand"]
    description: str = gene[1]['Description']
    if isinstance(description, str):
        for _ in "()`/,'":
            description: str = description.replace(_, "")
        description: str = description.replace(" ", "_")
    for fasta in SeqIO.parse((f"{busco_dir}/busco_sequences/"
                              f"single_copy_busco_sequences/"
                              f"{id}.fna"), "fasta"):
        seqs += f">{id}|{description}|{strand}|{fasta.id}\n"
        seqs += f"{fasta.seq}\n"

output = open("single_copy_busco_sequences.fasta", "w")
output.write(seqs)
output.close()
