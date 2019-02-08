
##########################################################
#
# CNV_Match.py
#
# Methods adopted from Connor Howington's VCF_Match.py
#
# Search all windows that have at least LOW_BOUND
# progenies with the value zero.
#
# Author: Fanli Zhou
# Date: 11/27/18
#
# Written using Python 3.6.5
#
##########################################################



import argparse, os, sys
from collections import defaultdict

def parse_args():

    parser = argparse.ArgumentParser(description='Takes a data file and finds windows of potential deletes. Output is printed to stdout.')
    parser.add_argument('file', help='an input data file', nargs='+')
    parser.add_argument('--low_bound', '-l', type=int, default = 5, help='only return windows that have at least LOW_BOUND zeros')

    args = parser.parse_args()
    for path in list(args.file):
        if not os.path.isfile(path):
            parser.error('File "{0}" cannot be found.'.format(path))

    
    return args

def find_windows(file_paths):

    # windows struction:
    # 1st dict key: chromosome -> 2nd dict key: window location -> count of zeros
    windows = defaultdict(lambda:defaultdict(int))
    for path in file_paths:
        file =  open(path, 'r')
        file.readline()

        for line in file:
            split_line = line.split()
            
            if not split_line[0].startswith('Pf3D7'):
                sys.stderr.write('Error in "{0}": Unexpected data format.\n'.format(path))
                sys.exit(1)
            loc = split_line[0].split('_')
            chrom = loc[0]+'_'+loc[1]+'_'+loc[2]
            pos = int(loc[-1])*300
            count = 0
            for num in split_line[1:]:
                if num == '0':
                    count += 1
			
            windows[chrom][pos] = count
    file.close()
            
    return windows


def find_deletions(file_paths, low_bound):

    windows = find_windows(file_paths)

    file = open('data.map', 'w')

    for chrom in sorted(windows.keys()):
        for pos in sorted(windows[chrom].keys()):
            if windows[chrom][pos] > low_bound:
                file.write(chrom+'\t'+str(pos)+'\n')

    file.close()

# Main flow

args = parse_args()

find_deletions(args.file, args.low_bound)      
