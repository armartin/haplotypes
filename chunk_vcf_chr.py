import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('--vcf')
parser.add_argument('--window_size', default = 10e6)
parser.add_argument('--window_overlap', default = 5e6)
parser.add_argument('--out')
args = parser.parse_args()

if args.vcf.endswith('gz'):
    vcf = gzip.open(args.vcf)
else:
    vcf = open(args.vcf)

chr_min=3e9
chr_max = 0
for line in vcf:
    line = line.strip()
    if line.startswith('##'):
        continue
    elif line.startswith('#CHROM'):
        vcf_header = line.split()
    else:
        chr_min = min(chr_min, int(line[1]))
        chr_max = max(chr_max, int(line[1]))

print [chr_min, chr_max]