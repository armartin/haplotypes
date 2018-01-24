import argparse
import gzip
import sys

from datetime import datetime


def current_time():
    return ' [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'


def open_file(filename):
    if filename.endswith('gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return my_file


def main(args):
    bim = open(args.bim)
    #  set up dicts to store necessary info about bim file
    print 'storing bim file info ' + current_time()
    pos_cm = {}
    pos_index = {}
    index_pos = {}
    index_count = {}
    counter = 0
    for line in bim:
        line = line.strip().split()
        index_pos[counter] = line[3]
        pos_index[line[3]] = counter
        pos_cm[line[3]] = line[2]
        index_count[counter] = 0
        counter += 1

    #  go through match file and count all overlaps
    match = open_file(args.match)
    print 'reading haplotype info' + current_time()
    i = 0
    for line in match:
        if i % 5000000 == 0:
            print 'line ' + str(i) + current_time()
        sys.stdout.flush()
        line = line.strip().split()
        start = pos_index[line[5]] #  get index of start position
        end = pos_index[line[6]]
        #  increment all positions spanned by a haplotype, e.g. range(x, y + 1)
        for interval in range(start, end+1):
            index_count[interval] += 1
        i += 1

    print 'writing output' + current_time()
    out = open(args.out, 'w')
    for pos in map(str, sorted(map(int, pos_cm.keys()))):
        out.write('\t'.join(map(str, [pos_index[pos], pos, pos_cm[pos], index_count[pos_index[pos]]])) + '\n')
    out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bim')
    parser.add_argument('--match')
    parser.add_argument('--out')

    args = parser.parse_args()
    main(args)