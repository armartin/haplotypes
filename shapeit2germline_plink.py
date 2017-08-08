"""
Converts phased shapeit files to germline plink input
"""

import argparse
import gzip
import makeMap
import itertools


def main(args):
    if args.haps.endswith('gz'):
        haps = gzip.open(args.haps)
    else:
        haps = open(args.haps)

    genmap = open(args.recomb_map)

    out_map = open(args.out + '.map', 'w')
    out_ped = open(args.out + '.ped', 'w')
    inds = open(args.sample)
    inds.readline()
    inds.readline()

    if args.keep is not None:
        keep = args.keep
        keep_inds = set()
        for line in keep:
            line = line.strip().split()
            keep_inds.add(line[1])

    haps_a = []
    haps_b = []
    for variant in makeMap.full_map(args.chr, genmap, haps, None, 'haps'):
        out_map.write('\t'.join(map(str, [variant['chr'], variant['rsid'], variant['gen_pos'], variant['phys_pos']])) + '\n')
        my_haps = map(lambda x: str(int(x) + 1), variant['haps'])
        haps_a.append(my_haps[0::2])
        haps_b.append(my_haps[1::2])
    out_map.close()

    ## this is a helper function to sort out haplotype merging orientation
    def roundrobin(*iterables):
        "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
        # Recipe credited to George Sakkis
        pending = len(iterables)
        nexts = itertools.cycle(iter(it).next for it in iterables)
        while pending:
            try:
                for next in nexts:
                    yield next()
            except StopIteration:
                pending -= 1
                nexts = itertools.cycle(itertools.islice(nexts, pending))

    flip_a = map(list, zip(*haps_a))
    flip_b = map(list, zip(*haps_b))
    combined_haps = zip(flip_a, flip_b)
    c = 0
    for ind in inds:
        ind = ind.strip().split()
        if ind[1] in keep_inds:
            hap_a = combined_haps[c][0]
            hap_b = combined_haps[c][1]
            final_haps = list(roundrobin(hap_a, hap_b))

            out_ped.write(ind[0] + ' ' + ind[1] + ' 0 0 0 -9 ')
            out_ped.write(' '.join(final_haps) + '\n')
        c += 1

    out_ped.close()

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
