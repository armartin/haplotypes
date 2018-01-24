"""
Converts vcf file to germline plink input

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
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  FR02_2  FR02_6
10      135656  rs10904561      G       T       .       PASS    .       GT      0|1     0|0
10      135708  rs7917054       A       G       .       PASS    .       GT      1|0     1|1

"""
import argparse
import gzip
import re
import makeMap

parser = argparse.ArgumentParser()
parser.add_argument('--vcf')
parser.add_argument('--recomb_map')
parser.add_argument('--out')
args = parser.parse_args()

if args.vcf.endswith('gz'):
    vcf = gzip.open(args.vcf)
else:
    vcf = open(args.vcf)

for line in vcf:
    line = line.strip()
    if line.startswith('##'):
        continue
    elif line.startswith('#CHROM'):
        vcf_header = line.split()

for variant in makeMap.full_map(chr, genmap, bim, map_bim):
    print variant
    #these are what are returned from full_map
    #[chr, rsid, str(proportion), phys_pos, a0, a1]