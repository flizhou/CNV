##############################################################
#
# find_deletion_windows.py
#
# Search at least LOW_BOUND consecutive windows with 10-40
# (number of progenies, by default) values below a cut-off as
# potential deletions.
#
# Author: Fanli Zhou
# Date: 11/21/18
#
# Written using Python 3.6.5
#
###############################################################

import argparse, io, sys, os
import collections

def parse_args():
    
    parser = argparse.ArgumentParser(description='Take a data file to find potential deletion windows.')
    parser.add_argument('data_file', help='an input data file', nargs = '+')
    parser.add_argument('--cut_off', '-c', type=float, required = False, default=0.01, help='Only condsiders values below CUT_OFF as potential deletions. Default = 0.01')
    parser.add_argument('--low_bound', '-l', type=int, required = False, default=3, help='Only returns results of at least LOW_BOUND consecutive windows. Default = 3')
    parser.add_argument('--interval', '-i', type=int, required = False, nargs = 2, help='Only considers windows with numbers of deletions within the interval. Default = [10, 40]. Example input: 10 40')
    
    args = parser.parse_args()
    paths = list(args.data_file)
    
    for path in paths:
        if not os.path.isfile(path):
            parser.error('File "{0}" cannot be found.'.format(path))

    if not args.interval:
        args.interval = [10, 40]

    return args
    
def read_data(file_paths):
    
    data = []
    info = []
    for path in file_paths:
        file = io.open(path)
        file.readline()

        for line in file:
            split_line = line.split()
            
            if not split_line[0].startswith('Pf3D7'):
                sys.stderr.write('Error in "{0}": Unexpected data format.\n'.format(path))
                sys.exit(1)
                
            if 'API' in split_line[0]:
                continue
            
            for n in range(1, len(split_line)):
                split_line[n] = float(split_line[n])

            info.append((split_line[0][:11], int(split_line[0][12:])*300))    
            data.append(split_line[1:])
            
        file.close()
        
    return (data, info)


def binary_insert(num, nums):

    l = 0
    r = len(nums)
    while l < r:
        m = l + (r-l)//2
        if num > nums[m]:
            l = m+1
        else:
            r = m
    nums.insert(r, num)

    
def find_cut_off(data, cut_off):

    win_num = len(data)                         
    pro_num = len(data[0])

    cut_off_num = int(win_num*cut_off)
    cut_off_values = []
    
    for pro in range(pro_num):

        # This keeps up to CUT_OFF_NUM windows sorted. After all the insertion, the last number in the list is the cut_off_value.
        below_cut = []
        
        for win in range(win_num):
            
            if len(below_cut) < cut_off_num:
                binary_insert(data[win][pro], below_cut)
            else:
                if data[win][pro] < below_cut[-1]:
                    below_cut.pop()
                    binary_insert(data[win][pro], below_cut)
        cut_off_values.append(below_cut[-1])
    return cut_off_values


def find_deletes(data, info, cut_off_values, low_bound, interval):

    win_num = len(data)                         
    pro_num = len(data[0])

    # This stores the interval of potential deletions as key and the number of progenies that contain the deletion as value.
    dic = collections.defaultdict(int)
    ret = {}
    for pro in range(pro_num):
        win = 0
        begin = end = -1
        cut_off = cut_off_values[pro]

        # In each progeny, find all the intervals with values no bigger than CUT_OFF.
        while win < win_num:                      
                                            
            if begin < 0 and data[win][pro] <= cut_off:
                begin = win
                
            elif begin >= 0:
                if data[win][pro] > cut_off:
                    end = win-1
                elif win == win_num-1:
                    end = win
                if end >= 0:
                    
                    # Find intervals that contain at least LOW_BOUND consecutive windows.
                    if end-begin >= low_bound-1:      
                        dic[(begin, end)] += 1
                    win = begin
                    begin = end = -1
            win += 1

    for key in dic:

        # Find intervals with numbers of deletions in the interval.
        if interval[0] < dic[key] < interval[1]:

            # Get the real chromosome location.
            loc_end = info[key[1]]         
            loc_beg = info[key[0]]

            # Combine all the intervals ending at the same place.
            if loc_end not in ret:          
                ret[loc_end]=(loc_beg, dic[key])
                
            elif loc_beg[1] < ret[loc_end][0][1]:
                ret[loc_end]=(loc_beg, dic[key])
                
    return ret


def write_in_file(data):
    
    file = open('deletion_windows.txt', 'w')
    
    file.write('deletion windows:\n')
    file.write('{0:11}  {1:8} {2:8} {3:8}\n'.format('Chromosome','begin', 'end', 'No(deletions)'))
    for end in sorted(data.keys()):
        file.write('{0:11}  {1:8} {2:8} {3:8}\n'.format(end[0], data[end][0][1], end[1]+300, data[end][1]))
                              
    file.close()
    
# Main flow
args = parse_args()
data, info = read_data(args.data_file)
cut_off_values = find_cut_off(data, args.cut_off)

deletes = find_deletes(data, info, cut_off_values, args.low_bound, args.interval)

write_in_file(deletes)
