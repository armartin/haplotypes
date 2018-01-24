"""
Converts HAPI-UR output to a phased VCF that can be used as input to BEAGLE for haplotype calling with e.g. RefinedIBD
"""

__author__ = 'armartin'

import argparse
from datetime import datetime
import time
import gzip

def open_gzipped(filename):
    """
    open potentially gzipped files
    """
    if filename.endswith('gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return(my_file)

def main(args):
    phgeno = open_gzipped(args.phgeno)
    phind = open_gzipped(args.phind)
    phsnp = open_gzipped(args.phsnp)
    
    out = gzip.open(args.out, 'w')
    out.write('##fileformat=VCFv4.2\n')
    out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
    
    ## write individual IDs
    for line in phind:
        line = line.strip().split()
        iid = line[0].split(':')[1].split('_')
        if line[0].endswith('_A'): #skip 1 per pair of haplotypes
            out.write('_'.join(iid[0:(len(iid)-1)]) + '\t')
    out.write('\n')
    
    ## write all SNPs
    for line in phsnp:
        ## write SNP info
        snp_line = line.strip().split()
        geno_line = phgeno.readline().strip()
        out.write('\t'.join([snp_line[1], snp_line[3], snp_line[0], snp_line[4], snp_line[5], '100', 'PASS', 'INFO', 'GT']))
        out.write('\t')
        
        ## write phased genotypes
        for i in zip(*(iter(range(len(geno_line))),) * 2): #get pairs at a time
            out.write(geno_line[i[0]] + '|' + geno_line[i[1]] + '\t')
        out.write('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--phgeno', required=True)
    parser.add_argument('--phind', required=True)
    parser.add_argument('--phsnp', required=True)
    
    parser.add_argument('--out', required=True)
    
    args = parser.parse_args()
    main(args)