#############################################
#
# find_deleted_genes.py
#
# Find genes in deletion windows
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
    parser.add_argument('window_file', help='The potential deletion windows file')
    parser.add_argument('genome_file', help='Gemone file from "http://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gff/data/"')

    args = parser.parse_args()

    for path in [args.window_file, args.genome_file]:
        if not os.path.isfile(path):
            parser.error('File "{0}" cannot be found.'.format(path))

    return args


def read_data(file_path):

    # deletions structure:
    # 1st dict key: chromosome -> 2nd dict key: window begin location -> 3rd dict key: window end location -> number of deletions
    deletions = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    file = io.open(file_path)
    file.readline()
    file.readline()

    for line in file:
        split_line = line.split()
        deletions[split_line[0]][int(split_line[1])][int(split_line[2])]=split_line[3]

    file.close()

    return deletions

     
class Gene_bank:

    def __init__(self, file_path):


        # gene_bank structure:
        # 1st dict key: chrom -> 2nd dict key: gene_begin location//10000 -> 3rd dict key: (gene_begin location, gene_end location) -> gene ID        
        self.gene_bank = defaultdict(lambda: defaultdict(dict))
        
        self.build_gene_bank(file_path)
        

    def build_gene_bank(self, file_path):
        
        file = io.open(file_path)
        
        for line in file:
            
            if line[0] == '#':
                continue

            if not line.startswith('Pf'):
                sys.stderr.write('Error in "{0}": Unexpected data format.\n'.format(file_path))
                sys.exit(1)
            
            split_line = line.split()
            if split_line[2] == 'gene':
                loc_chrom = split_line[0]
                loc_beg = split_line[3]
                loc_end = split_line[4]
                gene_id = split_line[-1]

                if not loc_beg or not loc_end:
                    continue
                
                loc_beg = int(loc_beg)
                loc_end = int(loc_end)
        
                code_beg = loc_beg//10000
                code_end = loc_end//10000     
                
                chrom = self.gene_bank[loc_chrom]
                chrom[code_beg][(loc_beg, loc_end)] = gene_id

                # If a gene spans the edges of two buckets, then it is added in both.
                if code_beg != code_end:                    

                    chrom[code_end][(loc_beg, loc_end)] = gene_id

        file.close()
        
        
    def search_genes(self, loc_chrom, loc_beg, loc_end):

        # genes structure:
        # genes key: gene location tuple(begin location, end location) -> gene ID
        genes = {}

        chrom = self.gene_bank[loc_chrom]
        code_beg = loc_beg//10000
        code_end = loc_end//10000
        
        if code_beg in chrom:
            for loc in chrom[code_beg]:
                
                if loc[0]>loc_end or loc[1]<loc_beg:
                    continue
                genes[loc] = chrom[code_beg][loc]
                
        if code_beg != code_end:
            if code_end in chrom:
                for loc in chrom[code_end]:
                    
                    if loc[0]>loc_end or loc[1]<loc_beg:
                        continue
                    genes[loc] = chrom[code_end][loc]
                    
        return genes

        
              
def find_deleted_genes(gene_bank, deletes):

    file = open('deleted_genes.txt', 'w')
    
    if deletes:
        file.write('deleted_genes:\n')
        file.write('{0:11}  {1:8} {2:8} {3}\n'.format('Chromosome', 'begin', 'end', 'gene ID'))
        
        for chrom in sorted(deletes.keys()):
            for loc_beg in sorted(deletes[chrom].keys()):
                for loc_end in sorted(deletes[chrom][loc_beg].keys()):
            
                    genes = gene_bank.search_genes(chrom, loc_beg, loc_end)

                    if genes:
                        
                        for loc in genes:
                             
                            file.write('{0:11}  {1:8} {2:8} {3}\n'.format(chrom, loc[0], loc[1], genes[loc]))
                    else:
                        file.write('No gene is found in this region\n\n')
    else:
        file.write("There's no deleted gene.")          



# Main flow
args = parse_args()
deletes = read_data(args.window_file)
gene_bank = Gene_bank(args.genome_file)
find_deleted_genes(gene_bank, deletes)
