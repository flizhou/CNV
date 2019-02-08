#############################################
#
# find_deleted_core_genes.py
#
# Find deleted core genes
#
# Author: Fanli Zhou
# Date: 11/21/18
#
# Written using Python 3.6.5
#
#############################################


import argparse, io, sys, os
from collections import defaultdict

def parse_args():

    parser = argparse.ArgumentParser(description='Find genes in the potential deletion windows.')
    parser.add_argument('gene_file', help='The deleted gene file')
    parser.add_argument('core_genome_file', help='Core_gemone file')

    args = parser.parse_args()

    for path in [args.gene_file, args.core_genome_file]:
        if not os.path.isfile(path):
            parser.error('File "{0}" cannot be found.'.format(path))

    return args


def read_data(file_path):

    # data struction:
    # 1st dict key: chromosome -> 2nd dict key: tuple(begin location, end location) -> gene ID or 'Core'
    data = defaultdict(lambda: defaultdict(tuple))

    file = io.open(file_path)

    for line in file:
        if line[0] != 'P':
            continue
        split_line = line.split()
        data[split_line[0]][(int(split_line[1]), int(split_line[2]))]=split_line[3]

    file.close()

    return data

def find_core_genes(core_genome, genes):
    
    file = open('deleted_core_genes.txt', 'w')

    if genes:
        file.write('deleted_core_genes:\n')
        file.write('{0:11}  {1:8} {2:8} {3}\n'.format('Chromosome', 'begin', 'end', 'gene ID'))

        for chrom in sorted(genes.keys()):
            for gene_loc in sorted(genes[chrom].keys()):
                if chrom in core_genome:
                    for core_loc in core_genome[chrom]:
                        if gene_loc[0] > core_loc[1] or gene_loc[1] < core_loc[0]:
                            continue
                        
                        file.write('{0:11}  {1:8} {2:8} {3}\n'.format(chrom, gene_loc[0], gene_loc[1], genes[chrom][gene_loc]))

    file.close()

     
# Main flow
args = parse_args()
genes = read_data(args.gene_file)
core_genome = read_data(args.core_genome_file)
find_core_genes(core_genome, genes)
