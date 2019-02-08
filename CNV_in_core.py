
######################################################
#
# CNV_in_core.py
#
# Methods adopted from Connor Howington's in-core.py
#
# Find deletion windows that are in core genome.
#
# Author: Fanli Zhou
# Date: 11/27/18
#
# Written using Python 3.6.5
#
######################################################

import argparse
from collections import defaultdict

def parse_args():
    
    parser = argparse.ArgumentParser(description='Determines whether a set of positions in a file are within core regions')

    parser.add_argument('core_file')
    parser.add_argument('file_path', help='an input data file', nargs='+')

    args = parser.parse_args()

    return args


def read_cores(path):
    
    # cores struction:
    # 1st dict key: chromosome -> [start pos, end pos, all information]
    cores = defaultdict(lambda: [])

    with open(path, 'r') as core_file:
        line = core_file.readline().rstrip()

        while line != '':

            split_line = line.split(' ')

            chrom = split_line[0]
            start = int(split_line[1])
            end = int(split_line[2])

            cores[chrom].append([start, end, line])
            line = core_file.readline()
        
    return cores


def in_core(file_path, cores):
    old_file = open(file_path, 'r')

    new_file = open(file_path+'.core', 'w', 1)

    line = old_file.readline().rstrip()

    while line != '':

        split_line = line.split('\t')

        chrom = split_line[0]
        pos = int(split_line[1])
        found = False

        if chrom in cores:
            for core in cores[chrom]:
                if core[0] <= pos <= core[1] or core[0] <= pos+300 <= core[1]:
                    found = True
                    new_file.write(line+'\t'+core[2].rstrip()+'\n')
                    break

        line = old_file.readline().rstrip()

    old_file.close()
    new_file.close()



# Main flow
args = parse_args()
cores = read_cores(args.core_file)

for path in args.file_path:
    in_core(path, cores)



