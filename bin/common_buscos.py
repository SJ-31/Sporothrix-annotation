import pandas as pd
from Bio import SeqIO

busco_dir = "-transcripts_merged.fasta_BUSCO/run_sordariomycetes_odb10"
common_genes: dict = {}
samples: tuple = ("S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10")

for sample in samples:
    directory = f"{sample}{busco_dir}"
    table = (f"{directory}/full_table.tsv")
    all_buscos = pd.read_csv(table, sep='\t', skiprows=2)
    complete = all_buscos[all_buscos['Status'] == 'Complete']
    for gene in complete.iterrows():
        id = gene[1]['# Busco id']
        description: str = gene[1]['Description']
        if isinstance(description, str):
            description: str = description.replace(" ", "_").replace("/", "")
        common_genes[id] = common_genes.get(id, [description, ''])
        for fasta in SeqIO.parse((f"{directory}/busco_sequences/"
                                  f"single_copy_busco_sequences/"
                                  f"{id}.fna"), "fasta"):
            common_genes[id][1] += f">{id}|{sample}|{description}|{fasta.id}\n"
            common_genes[id][1] += f"{fasta.seq}\n"

for gene_name, info in common_genes.items():
    write = open(f"{gene_name}-{info[0]}.all_samples.fasta", "w")
    write.write(f"{info[1]}")
    write.close()
