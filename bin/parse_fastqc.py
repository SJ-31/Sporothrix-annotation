import re
import subprocess

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

'''
Script for extracting modules from fastqc.txt output and doing visualizations
'''


def format(string: str) -> str:
    string = string.replace('>>', '').replace(' ', '_')
    string = re.sub('\t.*', '', string)
    return string.strip() + '.txt'


def get_num(string: str) -> int:
    return int(re.sub('[^0-9]+', '', string))


locs = subprocess.run(['cat fastqc_data | grep ">>" -n'],
                      shell=True,
                      stdout=subprocess.PIPE)
stdout = str(locs.stdout).split('\\n')[:-1]
modules = [get_num(num) for num in stdout if '_' not in num]
ends = [get_num(num) for num in stdout if 'END' in num]
with open('fastqc_data.txt') as f:
    to_parse = f.readlines()
    for strt, end in zip(modules, ends):
        print(format(to_parse[strt-1]))
        formatted = open(format(to_parse[strt-1]), 'w')
        formatted.write(''.join(to_parse[strt:end-1]))


# Basic statistics

# Per base sequence quality

# Per tile sequence quality

# Per base sequence content

# Per base N content

# Sequence length distribution

# Sequence duplication levels

# Overrepresented sequences

# Adapter content
