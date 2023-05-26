import sys

with open('genomic.gff') as f:
    lines = f.readlines()
new = []
for line in lines:
    if 'ID' in line:
        first_tab = line.find('\t')
        cut = f"{line[first_tab:]}"
        id = line[:line.find('\t')].replace('.1', '').replace('_', '')
        line = f'{id}{cut}'
        new.append(line)

sys.stdout.write(''.join(new))
