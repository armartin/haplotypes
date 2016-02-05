import argparse
import gzip
from scipy import stats
import numpy as np
import collections
import random
from bisect import bisect

## open gzipped and plain files
def open_file(filename):
    if filename.endswith('gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return my_file

## first step: perform t-test comparing in vs out of haplotype cluster
def true_test(dash, pheno_dict, all_inds):
    pheno_inds = set(pheno_dict.keys()).intersection(all_inds)
    in_clust_inds = set(dash[5:len(dash):2])
    
    in_clust_inds = list(pheno_inds.intersection(in_clust_inds).intersection(all_inds))
    not_in_clust_inds = list(pheno_inds.difference(in_clust_inds).intersection(all_inds))
    
    in_clust_phenos = [pheno_dict[ind] for ind in in_clust_inds]
    not_in_clust_phenos = [pheno_dict[ind] for ind in not_in_clust_inds]
    
    (my_t, my_p) = stats.ttest_ind(in_clust_phenos, not_in_clust_phenos)
    truth = {'in_clust': in_clust_inds, 'out_clust': not_in_clust_inds, 't': my_t, 'p': my_p}
    return truth

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
def perm_test(truth, ind_grid, pca_grid, pheno_dict, all_inds, times=100, p=[], t=[]):
    ## make sure there are greater than 2 individuals in and outside the cluster
    if len(truth['in_clust']) > 2:
        for i in range(times):
            matched_inds = []
            for ind in truth['in_clust']:
                matched = ind_grid[ind]
                matched_inds.append(random.sample(pca_grid[matched[0]][matched[1]], 1)[0]) #sample individual randomly
            unmatched_inds = list(set(truth['in_clust']).union(set(truth['out_clust'])).difference(set(matched_inds)))
            
            matched_phenos = [pheno_dict[ind] for ind in matched_inds]
            unmatched_phenos = [pheno_dict[ind] for ind in unmatched_inds]
            (my_t, my_p) = stats.ttest_ind(matched_phenos, unmatched_phenos)
            t.append(my_t)
            p.append(my_p)
        return({'t': sorted(t), 'p': sorted(p)})
    else:
        return({'t': 'NA', 'p': 'NA'})
    
## open all files
def main(args):
    dash = open_file(args.dash)
    pheno = open_file(args.pheno)
    #cum_ibd = open_file(args.cum_ibd)
    pca = open_file(args.pca)
    out = open(args.out, 'w')
    out.write('\t'.join(['cluster', 'cluster_start', 'cluster_end', 'max_start', 'min_end', 'cluster_midpoint', 'p_adj']) + '\n')
    
    all_inds = set()
    
    pheno_dict = {}
    for line in pheno:
        line = line.strip().split()
        try:
            pheno_dict[line[1]] = float(line[2]) #fix to make iid sometime
            all_inds.add(line[1])
        except ValueError:
            pass
    print '# of inds with phenos: '+ str(len(pheno_dict))
    
    pca_dict = {}
    pca_inds = set()
    pc1 = []
    pc2 = []
    for line in pca:
        line = line.strip().split()
        if line[1] in pheno_dict:
            pca_inds.add(line[1])
            pc1.append(float(line[2]))
            pc2.append(float(line[3]))
            pca_dict[line[1]] = map(float, line[2:len(line)])
    print '# of inds with PCs: ' + str(len(pca_dict))
    
    all_inds = list(set(pca_dict.keys()).intersection(set(pheno_dict.keys())))
    print '# of inds with PCs and phenos: ' + str(len(all_inds))
    
    ## make an evenly spaced pca grid for matching individuals, save individuals that fall in each grid
    
    pc1 = sorted(pc1)
    pc2 = sorted(pc2)
    pc1_grid = chunkIt(pc1, 10)
    pc1_bounds = [pc1_grid[i][0] for i in range(len(pc1_grid))]
    pc1_bounds.append(max(pc1))
    pc2_grid = chunkIt(pc2, 10)
    pc2_bounds = [pc2_grid[i][0] for i in range(len(pc2_grid))]
    pc2_bounds.append(max(pc2))
    pca_grid = collections.defaultdict(dict) #will need to change this to have perl's auto-vivification feature if we go deeper
    for i in range(len(pc1_grid)):
        for j in range(len(pc2_grid)):
            pca_grid[i][j] = set()
    
    print pc1_bounds
    print pc2_bounds
    
    ## store ind -> grid points
    ind_grid = {}
    for ind in all_inds:
        for i in range(len(pc1_bounds)-1):
            for j in range(len(pc2_bounds)-1):
                if pca_dict[ind][0] >= pc1_bounds[i] and pca_dict[ind][0] <= pc1_bounds[i+1] and pca_dict[ind][1] >= pc2_bounds[j] and pca_dict[ind][1] <= pc2_bounds[j+1]:
                    pca_grid[i][j].add(ind)
                    ind_grid[ind] = [i, j]
    
    #num inds in each grid point
    #for i in range(len(pc1_bounds)-1):
    #    for j in range(len(pc2_bounds)-1):
    #        print [i, j, len(pca_grid[i][j])]
    
    clust_dict = {}
    count = 1
    for line in dash:
        line = line.strip().split()
        clust_dict[line[0]] = line[1:5]
        truth = true_test(line, pheno_dict, all_inds)
        perm = perm_test(truth, ind_grid, pca_grid, pheno_dict, all_inds) #defaults to 100
        times = 100
        if perm['p'] != 'NA':
            p_adj = float(bisect(perm['p'], truth['p']))/len(perm['p'])
            print [count, p_adj]
            
            if p_adj <= 0.05:
                print 'more permutations'
                while p_adj < 5/times and times < 10001:
                    times = times * 10
                    perm = perm_test(truth, ind_grid, pca_grid, pheno_dict, times=times, p=perm['p'], t=perm['t']) #defaults to 100
                    p_adj = float(bisect(perm['p'], truth['p']))/len(perm['p'])
                    print [p_adj, times]
                    
            out.write('\t'.join(line[0:5]) + '\t' + str((int(line[1]) + int(line[2])) / 2) + '\t' + str(p_adj) + '\n')
        else:
            pass #maybe printing NA's is a good idea?
        count += 1
        if count % 1000 == 0:
            print [count, line[0:5]]
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dash', required=True)
    parser.add_argument('--cum_ibd')#, required=True)
    parser.add_argument('--pheno', help='this should be the joint maximal set of individuals phenotyped with haplotypes called', default='/homes/amartin/fr_broad/ldl_fr_engagex.pheno')
    parser.add_argument('--pca', required=True)    
    parser.add_argument('--out')#, required=True)
    
    args = parser.parse_args()
    main(args)