#!/bin/env Python
import sys
from Bio import SeqIO

types: dict = {'single_copy_busco_sequences': [],
               'multi_copy_busco_sequences': []}

look_in = sys.argv[1]
output_file_name = sys.argv[2]
find = set(sys.argv[3:])

busco_dir: str = f"{look_in}/run_sordariomycetes_odb10"
with open(f"{busco_dir}/full_table.tsv") as f:
    contents = f.read()
    records = contents[contents.rfind('#')+1:].splitlines()

for line in records:
    record = line.strip().split('\t')
    geneID = record[0]
    status = record[1]
    if status == 'Missing' and geneID in find:
        find.remove(geneID)
    elif status == 'Complete' and geneID in find:
        types['single_copy_busco_sequences'].append(geneID)
    elif status == 'Duplicated' and geneID in find:
        types['multi_copy_busco_sequences'].append(geneID)

result = ''
for t in types.keys():
    for gene in types[t]:
        for fasta in SeqIO.parse(f"{busco_dir}/busco_sequences/{t}/{gene}.fna",
                                 "fasta"):
            result += f">{gene}:{fasta.id}\n"
            result += f"{fasta.seq}\n"
final = open(f"{output_file_name}-found.fasta", "w")
final.write(result)
final.close()
