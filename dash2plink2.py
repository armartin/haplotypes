import argparse
import gzip
import numpy as np
from datetime import datetime

def open_file(filename):
    if filename.endswith('.gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return(my_file)

def current_time():
    return(' [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']')

def main(args):
    fam = open_file(args.fam)
    clst = open_file(args.clst)
    out_tped = open(args.out + '.tped', 'w')
    out_tfam = open(args.out + '.tfam', 'w')
    ind_order = []

    print 'Reading fam file for inds' + current_time()
    # write tfam inds and store in inds list
    for line in fam:
        line = line.strip().split()
        ind_order.append(line[1])
        out_tfam.write(' '.join(line) + '\n')
    out_tfam.close()

    ind_order = np.array(ind_order)

    print 'Reading clst file' + current_time()
    for line in clst:
        line = np.array(line.strip().split())
        out_tped.write(args.chr + ' ' + line[0] + ' 0 ' + str(int(line[1]) + int(line[2]) / 2) + ' ')
        ind_clusts = np.array(line[5:len(line):2])
        ind_chroms = np.array(line[6:len(line):2])
        haps = [1] * 2 * len(ind_order)
        count = [0] * 2 * len(ind_order)
        indices = np.in1d(ind_order, ind_clusts)

        indices = np.where(indices)
        diploid_indices = np.concatenate([indices[0]*2, indices[0]*2+1])
        for index in diploid_indices:
            if index % 2 == 0:
                if ind_order[index/2] + '.0' in ind_chroms:
                    haps[index] = 2
                    count[index] = 1
            else:
                if ind_order[index/2] + '.1' in ind_chroms:
                    haps[index] = 2
                    count[index] = 1
        out_tped.write(' '.join(map(str, haps)) + '\n')
    out_tped.close()

    print 'Done!' + current_time()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fam')
    parser.add_argument('--clst')
    parser.add_argument('--chr')
    parser.add_argument('--out')

    args = parser.parse_args()
    main(args)
