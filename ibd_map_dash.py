import argparse
import gzip
from scipy import stats

## open gzipped and plain files
def open_file(filename):
    if filename.endswith('gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return my_file

## first step: perform t-test comparing those with vs without haplotype
## compare pheno clust with pheno no-clust inds
## make sure there are >3 individuals in/out of cluster
def true_test(dash, pheno_dict):
    pheno_inds = set(pheno_dict.keys())
    in_clust = set(dash[5:len(dash):2])
    
    in_clust_pheno = pheno_inds.intersection(ind_clust)
    not_in_clust_pheno = pheno_inds.difference(in_clust_pheno)
    
    #print len(in_clust_pheno)
    #print len(not_in_clust_pheno)
    
    return(in_clust_pheno, not_in_clust_pheno)
    

## open all files
def main(args):
    dash = open_file(args.dash)
    pheno = open_file(args.pheno)
    #cum_ibd = open_file(args.cum_ibd)
    pca = open_file(args.pca)
    
    pheno_dict = {}
    for line in pheno:
        line = line.strip().split()
        pheno_dict[line[1]] = float(line[2])
    
    pca_dict = {}
    for line in pca:
        line = line.strip().split()
        pca_dict[line[1]] = map(float, line[6:len(line)])
    
    for line in dash:
        line = line.strip().split()
        (in_clust, not_in_clust) = true_test(line, pheno_dict)
        print len(in_clust)
        print len(not_in_clust)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dash', required=True)
    parser.add_argument('--cum_ibd')#, required=True)
    parser.add_argument('--pheno', help='this should be the joint maximal set of individuals phenotyped with haplotypes called', default='/homes/amartin/fr_broad/ldl_fr_engagex.pheno')
    parser.add_argument('--pca', required=True)    
    parser.add_argument('--out')#, required=True)
    
    args = parser.parse_args()
    main(args)