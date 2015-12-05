"""
Converts phased shapeit files to germline plink input

outputs
-------
map:
10	rs10904561	0.000099999997	135656
10	rs7917054	0.000099999997	135708
10	rs7089889	0.000099999997	178434

ped:
cas_dbs_db10_mix_jg_PSYC*PT-1BCE3 MMXII_iPSYCH_144731 0 0 2 2 1 1 2 1 1 1
cas_dbs_db10_mix_jg_PSYC*PT-1BCE4 MMXII_iPSYCH_144732 0 0 2 2 1 2 2 1 2 1
cas_dbs_db10_mix_jg_PSYC*PT-1BCE9 MMXII_iPSYCH_144737 0 0 2 2 1 2 2 1 2 1

input
-----
haps:
10 rs10904561 135656 G T 1 1 0
10 rs7917054 135708 A G 0 1 1 0
10 rs7089889 178434 T G 0 1 1 0

sample:
ID_1 ID_2 missing
0 0 0
egfext5401330 egfext5401330 0
egfext5401338 egfext5401338 0
egfext5401346 egfext5401346 0
"""

import argparse
import gzip
import makeMap
import itertools

parser = argparse.ArgumentParser()
parser.add_argument('--haps')
parser.add_argument('--sample')
parser.add_argument('--recomb_map')
parser.add_argument('--chr')
parser.add_argument('--out')
args = parser.parse_args()

if args.haps.endswith('gz'):
    haps = gzip.open(args.haps)
else:
    haps = open(args.haps)

genmap = open(args.recomb_map)

#for line in haps:
#    line = line.strip()
#    if line.startswith('##'):
#        continue
#    elif line.startswith('#CHROM'):
#        vcf_header = line.split()

out_map = open(args.out + '.map', 'w')
out_ped = open(args.out + '.ped', 'w')
inds = open(args.sample)
inds.readline()
inds.readline()
###NOTE: need to transpose whole haps file to write ped format

haps_a = []
haps_b = []
for variant in makeMap.full_map(args.chr, genmap, haps, None, 'haps'):
    #print variant
    out_map.write('\t'.join(map(str, [variant['chr'], variant['rsid'], variant['gen_pos'], variant['phys_pos']])) + '\n')
    my_haps = map(lambda x: str(int(x) + 1), variant['haps'])
    #print len(my_haps)
    #print type(my_haps)
    #print my_haps[0:10]
    haps_a.append(my_haps[0::2])
    haps_b.append(my_haps[1::2])
    #out_ped.write('\t'.join([variant['chr'], variant['rsid'], variant['gen_pos'], variant['phys_pos']]))
    #these are what are returned from full_map
    #[chr, rsid, str(proportion), phys_pos, a0, a1]

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
print len(flip_a)
#print len(flip_haps)
c = 0
for ind in inds:
    #print len(flip_a[c])
    hap_a = combined_haps[c][0]
    hap_b = combined_haps[c][1]
    final_haps = list(roundrobin(hap_a, hap_b))
    #print len(flip_haps[c])
    ind = ind.strip().split()
    out_ped.write(ind[0] + ' ' + ind[1] + ' 0 0 0 -9 ')
    out_ped.write(' '.join(final_haps) + '\n')
    c += 1
    
out_map.close()
out_ped.close()