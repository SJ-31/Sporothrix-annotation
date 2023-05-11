import argparse
import sys

import numpy as np

parser = argparse.ArgumentParser(description='Change maker control file'
                                 'params with'
                                 'command line arguments')
parser.add_argument('params', nargs='+', metavar='F_P',
                    help='F is the file the parameter is located in, '
                    'P the parameter (no spaces)\n'
                    'Ex: exe_RepeatMasker=/usr/bin/RepeatMasker/RepeatMasker '
                    'will change the location of the RepeatMasker executable'
                    'in the "maker_exe.ctl" file'
                    )
args = parser.parse_args()


def extr(string: str) -> tuple[str, str, str]:
    '''Returns a tuple of the '''
    filename: str = string[:string.find('_')]
    ori: str = string[string.find('_')+1:string.find('=')+1]
    change: str = string[string.find('_') + 1:]
    return filename, ori, change


def is_in(reference: list[str], query: list[str]) -> bool:
    '''Determine if a substring exists in a list of strings'''
    for sub in reference:
        if any(q in sub for q in query):
            return True
    return False


def change_param(filename: str, text: list[str], params: list[str]) -> str:
    changed: list[str] = []
    for line in text:
        move_on = False
        for index, param in enumerate(params):
            split_param = extr(param)
            if split_param[0] not in filename:
                continue
            if split_param[1] in line:
                params.pop(index)
                changed.append(f'{split_param[2]}\n')
                move_on = True
        if move_on:
            continue
        changed.append(line)
    return ''.join(changed)


files: dict[str, list[str]] = {}
for file in ['maker_evm.ctl', 'maker_exe.ctl',
             'maker_opts.ctl', 'maker_bopts.ctl']:
    with open(file) as f:
        files[file] = np.array(f.readlines())

all_params = sys.argv[1:]

param_queries = [extr(param)[1] for param in all_params]
for file, text in files.items():
    if is_in(text, param_queries):
        with open(file, 'w') as w:
            w.write(change_param(file, text, all_params))
