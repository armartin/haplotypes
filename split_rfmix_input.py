import argparse
import gzip
import numpy as np


#round to an arbitrary base
def myround(x, base=2):
    return int(base * round(float(x)/base))

#generate output data with admixed set split into groups
def main(args):
    if args.alleles.endswith('gz'):
        alleles = gzip.open(args.alleles)
    else:
        alleles = open(args.alleles)
    classes = open(args.classes).readline().strip().split()
    fam = open(args.fam)
    inds = []
    for line in fam:
        line = line.strip().split()
        inds.append(line)
    nsplits = int(args.nsplits)+1 #add 1, fence post problem for boundaries
    num_admixed_haps = classes.count('0')
    ref_indices = []
    admixed_indices = []
    for hap in range(len(classes)):
        if classes[hap] != '0':
            ref_indices.append(hap)
        else:
            admixed_indices.append(hap)
    
    print ref_indices
    #print admixed_indices
    
    boundaries = np.linspace(0,len(admixed_indices), nsplits)
    #require that an individual (2 haplotypes) are run together
    round_bound = [myround(bound) for bound in boundaries]
    
    print boundaries
    print round_bound
    
    prefix = args.out
    a_files = [open(prefix + '_run' + str(i) + '.alleles', 'w') for i in range(len(round_bound)-1)]
    c_files = [open(prefix + '_run' + str(i) + '.classes', 'w') for i in range(len(round_bound)-1)]
    f_files = [open(prefix + '_run' + str(i) + '.fam', 'w') for i in range(len(round_bound)-1)]
    for bound in range(len(round_bound)-1):
        [c_files[bound].write('0 ') for i in range(round_bound[bound], round_bound[bound+1])]
        [c_files[bound].write(classes[i] + ' ') for i in ref_indices]
        
    print ref_indices
    print len(round_bound)
    print len(inds)
    for bound in range(len(round_bound)-1):
        [f_files[bound].write(' '.join(inds[i]) + '\n') for i in range(round_bound[bound], round_bound[bound+1], 2)] #only do this every other
        [f_files[bound].write(' '.join(inds[i]) + '\n') for i in ref_indices[0::2]]
    
    for bound in range(len(round_bound)-1):
        [c_files[bound].close()]
        [f_files[bound].close()]
        
    for line in alleles:
        line = line.strip()
        for bound in range(len(round_bound)-1): # read through alleles once, write output for each run into separate files
            [a_files[bound].write(line[i]) for i in range(round_bound[bound], round_bound[bound+1])]
            [a_files[bound].write(line[i]) for i in ref_indices]
            a_files[bound].write('\n')
        

#parse arguments
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--alleles', required=True) #don't need snp_locations because these should all be the same
    parser.add_argument('--classes', required=True)
    parser.add_argument('--fam', required=True)
    parser.add_argument('--nsplits', help = 'number of ways the admixed individuals should be split (reference is constant)')
    parser.add_argument('--out')
    
    args = parser.parse_args()
    main(args)