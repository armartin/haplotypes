"""
Converts phased shapeit files to germline plink input

outputs:
-------
map:
10      rs10904561      0.000099999997  135656
10      rs7917054       0.000099999997  135708
10      rs7089889       0.000099999997  178434

ped - n x (6 + 2m), n=# inds, m=# markers:
cas_dbs_db10_mix_jg_PSYC*PT-1BCE3 MMXII_iPSYCH_144731 0 0 2 2 1 1 2 1 1 1
cas_dbs_db10_mix_jg_PSYC*PT-1BCE4 MMXII_iPSYCH_144732 0 0 2 2 1 2 2 1 2 1
cas_dbs_db10_mix_jg_PSYC*PT-1BCE9 MMXII_iPSYCH_144737 0 0 2 2 1 2 2 1 2 1

input
-----
haps.gz - m x (5 + 2n), n=# inds, m=#markers:
22 rs5746647 17057138 G T 1 1 1 1 1 1 0 1 1 1 1
22 rs5747988 17073066 A G 1 1 1 1 0 1 1 0 0 1 1
22 rs5747999 17075353 C A 1 1 1 1 1 1 1 0 0 1 1

sample:
ID_1 ID_2 missing father mother sex plink_pheno
0 0 0 D D D B
1 1 0 0 0 2 -9
2 2 0 0 0 1 -9
3 3 0 0 0 2 -9
"""

import argparse
import gzip
import makeMap
from datetime import datetime
import itertools


def current_time():
    return ' [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'


def open_file(filename):
    if filename.endswith('.gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return my_file


def roundrobin(*iterables):
    """roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    this is a helper function to sort out haplotype merging orientation
    Recipe credited to George Sakkis"""
    pending = len(iterables)
    nexts = itertools.cycle(iter(it).next for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = itertools.cycle(itertools.islice(nexts, pending))

def main(args):
    haps = open_file(args.haps)
    genmap = open_file(args.recomb_map)

    out_map = open(args.out + '.map', 'w')
    out_ped = open(args.out + '.ped', 'w')
    inds = open(args.sample)
    inds.readline()
    inds.readline()

    if args.keep is not None:
        keep = open(args.keep)
        keep_inds = set()
        for line in keep:
            keep_inds.add(line.strip().split()[1])

    haps_a = []
    haps_b = []
    print 'Writing map' + current_time()
    for variant in makeMap.full_map(args.chr, genmap, haps, None, 'haps'):
        out_map.write('\t'.join(map(str, [variant['chr'], variant['rsid'], variant['gen_pos'], variant['phys_pos']])) + '\n')
        my_haps = map(lambda x: str(int(x) + 1), variant['haps'])
        haps_a.append(my_haps[0::2])
        haps_b.append(my_haps[1::2])
    out_map.close()

    print 'Flipping haps' + current_time()
    # flipping dimension (sort of like transpose)
    flip_a = map(list, zip(*haps_a))
    flip_b = map(list, zip(*haps_b))
    combined_haps = zip(flip_a, flip_b)

    c = 0
    inds = open(args.sample)
    inds.readline()
    inds.readline()
    for ind in inds:
        ind = ind.strip().split()
        if args.keep is None or args.keep is not None and ind[1] in keep_inds:
            hap_a = combined_haps[c][0]
            hap_b = combined_haps[c][1]
            final_haps = list(roundrobin(hap_a, hap_b))

            out_ped.write(ind[0] + ' ' + ind[1] + ' 0 0 0 -9 ')
            out_ped.write(' '.join(final_haps) + '\n')
        c += 1

    out_ped.close()
    print 'Done!' + current_time()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--haps')
    parser.add_argument('--sample')
    parser.add_argument('--keep')
    parser.add_argument('--recomb_map')
    parser.add_argument('--chr')
    parser.add_argument('--out')
    args = parser.parse_args()

    main(args)
