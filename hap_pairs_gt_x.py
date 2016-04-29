"""
Counts number of haplotypes per pair greater than length X
"""

import argparse
import gzip
import numpy as np

def main(args):
    pair_count = 0
    cum_pairs = gzip.open(args.cum_pairs)
    bins = np.arange(1,10.1,0.1)
    num_pairs = 0
    hist_bins = {}
    for b in bins:
        hist_bins[b] = 0
    for line in cum_pairs:
        num_pairs += 1
        line = line.strip().split()
        hap_sizes = sorted(map(float, line[3]))
        for b in hist_bins:
            hist_bins[b] += sum(i > b for i in hap_sizes)
        print num_pairs
        print hist_bins
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cum_pairs', default='/home/unix/armartin/finrisk/haplotypes/b37_gaps.txt')
    parser.add_argument('--out')
    
    args = parser.parse_args()
    main(args)
