import argparse
import gzip
import numpy as np

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
        line = line.split()
        chr_min = min(chr_min, int(line[1]))
        chr_max = max(chr_max, int(line[1]))
        chr = line[0]

print [chr_min, chr_max]

windows1 = np.arange(chr_min, chr_max, args.window_size)
windows2 = np.arange(chr_min+args.window_overlap, chr_max+args.window_overlap, args.window_size)
if windows2[-1] < chr_max:
    windows2[-1] = chr_max #make sure all SNPs make it in at least one window

print windows1
print windows2

out = open(args.out, 'w')
for i in range(len(windows1)-1):
    out.write(chr + ':' + str(int(windows1[i])) + '-' + str(int(windows1[i+1])) + '\n')
    out.write(chr + ':' + str(int(windows2[i])) + '-' + str(int(windows2[i+1])) + '\n')
out.close()