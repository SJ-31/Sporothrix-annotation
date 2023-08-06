#!/usr/bin/env python
# Script writing ragout rcp file

import sys

params = sys.argv[1:]
query = params[0]
references = params[1:]

recipe = []
formatted = query[:query.find("_")]
recipe.append(f'.target = {formatted}\n')
ref_list = '.references = '
ref_list += ','.join([r.replace(".fasta", "") for r in references])
recipe.extend([ref_list, '\n\n', f'{formatted}.fasta = ./{query}\n'])
for ref in references:
    recipe.append(f'{ref} = ./{ref}\n')
with open('recipe.rcp', 'w') as r:
    r.write(''.join(recipe))
