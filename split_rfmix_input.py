import argparse
import gzip
import numpy as np


#round to an arbitrary base
def myround(x, base=2):
    return int(base * round(float(x)/base))

#generate output data with admixed set split into groups
def main(args):
    alleles = open(args.alleles)
    classes = open(args.classes).readline().strip().split()
    print classes
    nsplits = args.nsplits #note: require that an individual (2 haplotypes) are run together
    num_admixed_haps = classes.count('0')
    ref_indices = []
    admixed_indices = []
    for hap in classes:
        if classes[hap] != '0':
            ref_indices.append(hap)
        else:
            admixed_indices.append(hap)
    
    print ref_indices
    print admixed_indices
    
    boundaries = np.arange(0,len(admixed_indices), step=nsplits)
    round_bound = [myround(bound) for bound in boundaries]
    
    print boundaries
    print round_bound

#parse arguments
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--alleles', required=True) #don't need snp_locations because these should all be the same
    parser.add_argument('--classes', required=True)
    parser.add_argument('--nsplits', help = 'number of ways the admixed individuals should be split (reference is constant)')
    parser.add_argument('--out')
    
    args = parser.parse_args()
    main(args)