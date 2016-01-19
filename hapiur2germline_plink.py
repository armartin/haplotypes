"""
Converts phased hapi-ur files to germline plink input

outputs
-------
map:
10	rs10904561	0.000099999997	135656
10	rs7917054	0.000099999997	135708
10	rs7089889	0.000099999997	178434

ped - n x (6 + 2m), n=# inds, m=# markers:
cas_dbs_db10_mix_jg_PSYC*PT-1BCE3 MMXII_iPSYCH_144731 0 0 2 2 1 1 2 1 1 1
cas_dbs_db10_mix_jg_PSYC*PT-1BCE4 MMXII_iPSYCH_144732 0 0 2 2 1 2 2 1 2 1
cas_dbs_db10_mix_jg_PSYC*PT-1BCE9 MMXII_iPSYCH_144737 0 0 2 2 1 2 2 1 2 1

input
-----
phsnp:
    rs8142331  22        2.125947237015        16494187 A C
    rs5746647  22        2.340497732162        17057138 G T
   rs16980739  22        2.344379186630        17058616 T C

fam:
cas_pts_cogy_aam_ld_Il5M*0 cas_pts_cogy_aam_ld_Il5M*0:122020511 0 0 2 2
cas_pts_cogy_aam_ld_Il5M*0 cas_pts_cogy_aam_ld_Il5M*0:122031238 0 0 2 2
cas_pts_cogy_aam_ld_Il5M*0 cas_pts_cogy_aam_ld_Il5M*0:122041977 0 0 2 2

phgeno - m x 2n:
00010000000000100000
11100001111000001111
01000000000001000001
"""

import argparse
import gzip
import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--phgeno')
parser.add_argument('--phsnp')
parser.add_argument('--fam')
parser.add_argument('--out')
args = parser.parse_args()

## open file, regardless of gzipped status
def open_file(filename):
    if filename.endswith('gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return my_file

phgeno = open_file(args.phgeno)
phsnp = open_file(args.phsnp)

out_map = open(args.out + '.map', 'w')
out_ped = open(args.out + '.ped', 'w')
inds = open(args.fam)

for line in phsnp:
    line = line.strip().split()
    out_map.write('\t'.join([line[1], line[0], line[2], line[3]]) + '\n')
out_map.close()

haps_a = []
haps_b = []
for line in phgeno:
    line = line.strip().split()
    my_haps = map(lambda x: str(int(x) + 1), line)
    haps_a.append(my_haps[0::2])
    haps_b.append(my_haps[1::2])

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
print len(flip_a)
print len(flip_a[0])
combined_haps = zip(flip_a[0], flip_b[0])
c = 0
for ind in inds:
    ind = ind.strip().split()
    #print combined_haps
    print combined_haps[c]
    hap_a = combined_haps[c][0]
    hap_b = combined_haps[c][1]
    final_haps = list(roundrobin(hap_a, hap_b))
    ind = ind.strip().split()
    out_ped.write(' '.join(ind) + '\t')
    out_ped.write(' '.join(final_haps) + '\n')
    c += 1
out_ped.close()
