import sys
# Remove unwanted characters for fasta for use with Augustus
with open('GCF_000961545.1_S_schenckii_v1_genomic.fna') as f:
    lines = f.readlines()
new = []
for line in lines:
    if line.startswith(">"):
        cleaned = line[:line.find(' ')]
        cleaned = cleaned.replace('.1', '').replace('_', '')
        new.append(f'{cleaned}\n')
        continue
    new.append(line)

sys.stdout.write(''.join(new))
