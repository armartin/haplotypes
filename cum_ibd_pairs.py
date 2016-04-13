import argparse
import gzip
import collections
from datetime import datetime

def main(args):
    if args.match.endswith('gz'):
        match = gzip.open(args.match)
    else:
        match = open(args.match)
    
    ## first, make sure all pairs are represented
    sample = open(args.sample)
    sample.readline()
    sample.readline()
    samples = []
    for line in sample:
        line = line.strip().split()
        samples.append(line[1])
    
    cum_ibd = {}
    for ind1 in samples:
        for ind2 in samples:
            if ind1 != ind2:
                ind_pairs = tuple(sorted((ind1, ind2)))
                cum_ibd[ind_pairs] = 0
    
    ## next, count cumulative IBD sharing across pairs
    for line in match:
        line = line.strip().split()
        ind1 = line[1].split('.')[0]
        ind2 = line[3].split('.')[0]
        if ind1 != ind2:
            ind_pairs = tuple(sorted((ind1, ind2)))
            cum_ibd[ind_pairs] += float(line[10])
    
    ## test
    all_pairs = cum_ibd.keys()
    print all_pairs[0]
    print cum_ibd[all_pairs[0]]
    print
    print all_pairs[1]
    print cum_ibd[all_pairs[1]]
    print


if __name__ == '__main__':    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--match')
    parser.add_argument('--sample')
    parser.add_argument('--out')
    
    args = parser.parse_args()
    main(args)