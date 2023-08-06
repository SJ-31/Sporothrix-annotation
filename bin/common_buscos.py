#!/bin/env python
import os
import pandas as pd
lineage = "sordariomycetes_odb10"
# Run this in the directory containing the busco folders
look_in = list(filter(lambda d: os.path.isdir(d), os.listdir(".")))
gene_dict: dict = {}
for directory in look_in:
    path = f"{directory}/run_{lineage}/full_table.tsv"
    busco_table = pd.read_csv(path, sep="\t", skiprows=2)
    df = busco_table[busco_table['Status'] == 'Complete']
    gene_dict[directory] = set(df['# Busco id'])

common_genes: set = set.intersection(*list(gene_dict.values()))
with open('common_genes.txt', 'w') as g:
    g.write("".join([f"{gene}\n" for gene in common_genes]))
