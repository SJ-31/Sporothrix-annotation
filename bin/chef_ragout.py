#!/usr/bin/env python

import sys

params = sys.argv[1:]
query = params[0]
references = params[1:]

recipe = []
recipe.append(f'.target = {query[:query.find("_")]}\n')
ref_list = '.references = '
ref_list += ','.join([r.replace(".fasta", "") for r in references])
recipe.extend([ref_list, '\n\n', f'{query} = ./{query}\n'])
for ref in references:
    recipe.append(f'{ref} = ./{ref}\n')
with open('recipe.rcp', 'w') as r:
    r.write(''.join(recipe))
