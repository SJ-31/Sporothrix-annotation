#!/bin/env python
import pandas as pd

look_in = [f"S0{i}" for i in range(2, 10)] + ["S10"]
gene_dict: dict = {}
for directory in look_in:
    path = f"{directory}_scaffolds_BUSCO/run_sordariomycetes_odb10/full_table.tsv"
    busco_table = pd.read_csv(path, sep="\t", skiprows=2)
    df = busco_table[busco_table['Status'] == 'Complete']
    gene_dict[directory] = set(df['# Busco id'])

common_genes: set = set.intersection(*list(gene_dict.values()))
with open('common_genes.txt', 'w') as g:
    g.write("".join([f"{gene}\n" for gene in common_genes]))
