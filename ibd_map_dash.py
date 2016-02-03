import argparse
import gzip
from scipy import stats
import numpy as np
import collections

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
    in_clust_inds = set(dash[5:len(dash):2])
    
    in_clust_inds = list(pheno_inds.intersection(in_clust_inds))
    not_in_clust_inds = list(pheno_inds.difference(in_clust_inds))
    
    in_clust_phenos = [pheno_dict[ind] for ind in in_clust_inds]
    not_in_clust_phenos = [pheno_dict[ind] for ind in not_in_clust_inds]
    
    if len(in_clust_phenos) > 2 and len(not_in_clust_phenos) > 2:
        print 'in test scenario'
        (my_t, my_p) = stats.ttest_ind(in_clust_phenos, not_in_clust_phenos)
        return {'in_clust': in_clust_inds, 'out_clust': not_in_clust_inds, 't': my_t, 'p': my_p}
        #print len(in_clust_pheno)
    #print len(not_in_clust_pheno)
    
    return(in_clust_phenos, not_in_clust_phenos)

## function for splitting into evenly sized grid
def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg  
    return out

## second step: perform t-test comparing individuals matched to those with vs without haplotype
def perm_test(truth, times=100):
    pass

## open all files
def main(args):
    dash = open_file(args.dash)
    pheno = open_file(args.pheno)
    #cum_ibd = open_file(args.cum_ibd)
    pca = open_file(args.pca)
    
    pheno_dict = {}
    for line in pheno:
        line = line.strip().split()
        try:
            pheno_dict[line[0]] = float(line[2]) #fix to make iid sometime
        except ValueError:
            pass
    print len(pheno_dict)
    
    pca_dict = {}
    pc1 = []
    pc2 = []
    for line in pca:
        line = line.strip().split()
        if line[1] in pheno_dict:
            pc1.append(float(line[6]))
            pc2.append(float(line[7]))
            pca_dict[line[1]] = map(float, line[6:len(line)])
    
    ## make an evenly spaced pca grid for matching individuals
    print len(pca_dict)
    print pca_dict['FR97_2347']
    print
    pc1 = sorted(pc1)
    pc2 = sorted(pc2)
    pc1_grid = chunkIt(pc1, 10)
    pc1_bounds = [pc1_grid[i][0] for i in range(len(pc1_grid))]
    pc2_grid = chunkIt(pc2, 10)
    pc2_bounds = [pc2_grid[i][0] for i in range(len(pc2_grid))]
    pca_grid = collections.defaultdict(dict) #will need to change this to have perl's auto-vivification feature if we go deeper
    for i in range(len(pc1_grid)):
        for j in range(len(pc2_grid)):
            pca_grid[i][j] = set()
    
    print pc1_bounds
    print pc2_bounds
    
    #print pca_grid
    for ind in pheno_dict.keys():
        for i in range(len(pc1_bounds)-1):
            for j in range(len(pc2_bounds)-1):
                if pca_dict[ind][0] >= pc1_bounds[i] and pca_dict[ind][0] < pc1_bounds[i+1] and pca_dict[ind][1] >= pc2_bounds[j] and pca_dict[ind][1] < pc2_bounds[j+1]:
                    pca_grid[i][j].add(ind)
                #else:
                #    print [ind, pca_dict[ind], i, j]
                #except KeyError:
                #    print [ind, pca_dict[ind], i, j]
    for i in range(len(pc1_bounds)-1):
        for j in range(len(pc2_bounds)-1):
            print len(pca_grid[i][j])
    
    
    clust_dict = {}
    for line in dash:
        line = line.strip().split()
        clust_dict[line[0]] = line[1:5]
        truth = true_test(line, pheno_dict)
        perm = perm_test(truth, times) #defaults to 100
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dash', required=True)
    parser.add_argument('--cum_ibd')#, required=True)
    parser.add_argument('--pheno', help='this should be the joint maximal set of individuals phenotyped with haplotypes called', default='/homes/amartin/fr_broad/ldl_fr_engagex.pheno')
    parser.add_argument('--pca', required=True)    
    parser.add_argument('--out')#, required=True)
    
    args = parser.parse_args()
    main(args)