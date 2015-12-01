import argparse
import re
import gzip
import collections
from datetime import datetime
import time
import numpy as np
from scipy import stats
import random

parser = argparse.ArgumentParser()
parser.add_argument('--ibd', default='/fs/projects/popgen_genotypes/haplotypes/ibdseq/ibd/Engagex_Egfext_Migraine_FinnishCardio_FTC_chr22.ibdseq.ibd2.gz')
parser.add_argument('--cum_ibd', default='/homes/amartin/popgen_genotypes/haplotypes/ibdseq/ibd/Engagex_Egfext_Migraine_FinnishCardio_FTC.sorted.cum_pairs.ibd.gz')
parser.add_argument('--pos', default='/homes/amartin/popgen_genotypes/IlluminaCoreExome/Engagex_Egfext_Migraine_FinnishCardio_FTC_projects/shapeit2_haplotypes/Engagex_Egfext_Migraine_FinnishCardio_FTC_chr22.haps.gz')
parser.add_argument('--pheno', help='this should be the joint maximal set of individuals phenotyped with haplotypes called', default='/homes/amartin/fr_broad/ldl_fr_engagex.pheno')
#parser.add_argument('--sample', default='/homes/amartin/popgen_genotypes/IlluminaCoreExome/all_birth_wgs.csv')
parser.add_argument('--pca', default='/homes/amartin/popgen_genotypes/haplotypes/pca/Palotie_Engagex_geno05_maf05_ld50.pcs')

parser.add_argument('--out', default='/fs/projects/popgen_genotypes/haplotypes/ibdseq/ibd/Engagex_Egfext_Migraine_FinnishCardio_FTC_chr22.hapassoc')
args = parser.parse_args()
""" parameters to consider:
    minimum permutation threshold
    case/control or continuous
    how to learn SNP resolution from data given IBD density differences (e.g. HLA vs centromere)
    haplotype window to consider (chunks)
    use all inds in pheno file unless a keep samples file is specified
        keep samples
    column indices for ibd, cum_ibd, pos
    ibd LOD or quality threshold
    maybe allow multiple phenotypes?
    different test for different phenotype distributions (default: rank sum)
    consider using recombination map instead of physical positions
    only checking for excess of IBD? anything about haplotype structure?
"""

pc1_bin = np.arange(-5, 6)
pc2_bin = np.arange(-3.5, 4, 0.5)
pca_bin = collections.defaultdict(dict)

## identifies bins
def which_bin(perm_bins, current_value):
    if current_value > perm_bins[-1]:
        current_bin = len(perm_bins)
    elif current_value < perm_bins[0]:
        current_bin = 0
    else:
        current_bin = np.nonzero(perm_bins >= current_value)[0][0]
    return current_bin

#finds bins for each pc
ind_pca_bin = {}
for pc1 in range(len(pc1_bin)+1):
    for pc2 in range(len(pc2_bin)+1):
        pca_bin[pc1][pc2] = set()


#add inds to pc dict, store bins
pca = open(args.pca)
for line in pca:
    line = line.strip().split()
    ind_pc1 = float(line[6])
    ind_pc2 = float(line[7])
    bin1 = which_bin(pc1_bin,ind_pc1)
    bin2 = which_bin(pc2_bin,ind_pc2)
    #bin1 = pc1_bin[(ind_pc1>=pc1_bin) & (ind_pc1<=pc1_bin)]
    #bin2 = pc2_bin[(ind_pc2>=pc2_bin) & (ind_pc2<=pc2_bin)]
    pca_bin[bin1][bin2].add(line[1])
    ind_pca_bin[line[1]] = [bin1, bin2]

## function to open files and read in headers
def read_header(filename, is_gzipped=False, is_csv=False):
    if is_gzipped == True:
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    if is_csv:
        line = my_file.readline().strip().split(',')
    else:
        line = my_file.readline().strip().split()
    header_dict = {}
    for i in range(len(line)):
        header_dict[line[i]] = i
    return(my_file, header_dict)


## read in cumulative shared IBD across pairs
#(cum_ibd, cum_ibd_header) = read_header(args.cum_ibd, True)
perm_bins = np.arange(0, 150, 5) #define step size throughout code

#cum_ibd_pairs = {} # (ind1, ind2) -> cum ibd
#these are backwards shitty names. rename
#cum_ibd_bin = {} # bin -> (ind1, ind2)
#bin_cum_ibd = {} # (ind1, ind2) -> bin, where IDs are finrisk_ids
#print 'Reading IBD pairs and making permutation bins [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
#this currently takes ~4 minutes to read chr22 with 11,639 finns
#for line in cum_ibd:
#    line = line.strip().split()
#    ind1 = line[cum_ibd_header['FR1']]
#    ind2 = line[cum_ibd_header['FR2']]
#    inds = tuple(sorted([ind1, ind2]))
#    cum_ibd_pair = float(line[cum_ibd_header['cum_ibd']])/1e6
#    cum_ibd_pairs[inds] = cum_ibd_pair
#    #find bin
#    current_bin = which_bin(perm_bins, cum_ibd_pair)
#    #put inds in dict of bins
#    if current_bin in cum_ibd_bin:
#        cum_ibd_bin[current_bin].append(inds)
#    else:
#        cum_ibd_bin[current_bin] = [inds]
#    bin_cum_ibd[inds] = current_bin
    #if ibd_bin 
#TO DO: afterwards, cum_ibd_bin values into np.array for every key for faster sampling

### read in samples and make dict map for IDs
#print 'Reading sample ID mapping [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
#(samples, sample_header) = read_header(args.sample, False, True)
#sample_dict = {}
#all_inds = set()
#for line in samples:
#    line = line.strip().split(',')
#    ind = line[sample_header['Finrisk_ID']]
#    if line[sample_header['Project.name']] == 'Engagex_Egfext_Migraine_FinnishCardio_FTC':
#        all_inds.add(ind)
#    sample_dict[line[sample_header['ID2']]] = ind

## read in phenotypes. TO DO: write an option for case/control vs quantitative trait
print 'Reading all phenotypes [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
pheno = open(args.pheno)
pheno_dict = {}
pheno_inds = set()
for line in pheno:
    line = line.strip().split()
    try:
        #pheno_dict[sample_dict[line[1]]] = float(line[2]) #allow this to be int in c/c
        pheno_dict[line[0]] = float(line[2]) #allow this to be int in c/c
    except ValueError:
        #pheno_dict[sample_dict[line[1]]] = None #this might be a pain
        pheno_dict[line[0]] = None #this might be a pain
    #pheno_inds.add(sample_dict[line[1]])
    pheno_inds.add(line[0])
inds = sorted(pheno_dict.keys())
phenos = sorted(pheno_dict.values())
#print inds[0:10]
#print phenos[len(phenos)-10:len(phenos)]
ind_dict = {}
for ind in range(len(inds)):
    ind_dict[inds[ind]] = ind

## open positions, set up position variables
print 'Reading positions [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
pos = gzip.open(args.pos)
pos_ibd_phenos = {} #keep track of pheno pairs per site with ibd
pos_ibd_inds = {} #keep track of ind pairs per site with ibd
pos_added = {}
for line in pos:
    line = line.strip().split()
    site = int(line[2])
    pos_ibd_inds[site] = set()
    pos_ibd_phenos[site] = []
    pos_added[site] = set()
pos_array = np.array(sorted(pos_added.keys()))

## read in ibd segments shared across pairs, storing individuals and phenotypes for every position
print 'Reading IBD segments are storing positional IBD spanning [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
ibd = gzip.open(args.ibd)
#ibd_starts = collections.defaultdict(dict)
#ibd_ends = collections.defaultdict(dict)
print 'Reading IBD pairs [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
c=0
for line in ibd:
    if c % 10000 == 0:
        print str(c) + ': [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    c += 1
    line = line.strip().split()
    ind1 = line[0]
    ind2 = line[2]
    inds = tuple(sorted([ind1, ind2]))
    start = int(line[5])
    end = int(line[6])
    #if ind1 not in sample_dict or ind2 not in sample_dict:
    #    continue #not sure why this isn't working
    all_pos = pos_array[(pos_array>=start) & (pos_array<=end)] #this is probably most computationally expensive
    for current_pos in all_pos:
        #pheno1 = pheno_dict[sample_dict[ind1]]
        #pheno2 = pheno_dict[sample_dict[ind2]]
        if ind1 in pheno_dict and ind2 in pheno_dict:
            if pheno_dict[ind1] is not None and pheno_dict[ind2] is not None: #subset inds and then wouldn't have to deal with this
                pos_ibd_inds[current_pos].add(inds)
                pheno1 = pheno_dict[ind1]
                pheno2 = pheno_dict[ind2]
                
                pos_ibd_phenos[current_pos].append(pheno1)
                pos_ibd_phenos[current_pos].append(pheno2)
        #except KeyError:
            #continue
        
        #deprecate
        #try: # a bit more computationally expensive but also more powerful for continuous than c/c traits (more permutations)
        #    #TO DO: add entire tuple
        #    if sample_dict[ind1] not in pos_added[current_pos]:
        #        pos_ibd_phenos[current_pos].append(pheno_dict[sample_dict[ind1]])
        #    if sample_dict[ind2] not in pos_added[current_pos]:
        #        pos_ibd_phenos[current_pos].append(pheno_dict[sample_dict[ind2]])
        #    pos_added[current_pos].add(sample_dict[ind1]) #lengths at a given position are different for pos_added and pos_ibd_phenos
        #    pos_added[current_pos].add(sample_dict[ind2]) #why?
        #    ####NOTE! Check if I've already added individual to pos_added
    #if ind2 in ibd_starts[ind1]:
    #    ibd_starts[ind1][ind2].append(start)
    #    ibd_ends[ind1][ind2].append(end)
    #else:
    #    ibd_starts[ind1][ind2] = [start]
    #    ibd_ends[ind1][ind2] = [end]

## run association and permutations for each site
#cum_ibd_pairs = {} # (ind1, ind2) -> cum ibd
##these are backwards shitty names. rename
#cum_ibd_bin = {} # bin -> (ind1, ind2)
#bin_cum_ibd = {} # (ind1, ind2) -> bin

#instead of permuting on pairs of cumulative IBD, permute on PCA grid
def run_perms(permutations, mean_perms):
    for perm in range(permutations):
        perm_inds = []
        perm_pheno = []
        for my_bin in range(len(true_bins)):
            current_bin = true_bins[my_bin]
            perm_ind_pair = random.choice(cum_ibd_bin[current_bin])
            while pheno_dict[perm_ind_pair[0]] is None or pheno_dict[perm_ind_pair[1]] is None:
                perm_ind_pair = random.choice(cum_ibd_bin[current_bin])
            #need to make sure we match the same number as in true_phenos, might need to this several times
            #while perm_ind_pair[0] not in sample_dict or perm_ind_pair[1] not in sample_dict:
            #perm_pheno.append(pheno_dict[sample_dict[perm_ind_pair[0]]])
            #perm_pheno.append(pheno_dict[sample_dict[perm_ind_pair[1]]])
            perm_pheno.append(pheno_dict[perm_ind_pair[0]])
            perm_pheno.append(pheno_dict[perm_ind_pair[1]])
        mean_perms.append(np.mean(perm_pheno))
    return mean_perms

def run_perms_pca(permutations, mean_perms):
    for perm in range(permutations):
        perm_inds = []
        perm_pheno = []
        for my_bin in range(len(true_bins)):
            current_bin = true_bins[my_bin]
            perm_ind = random.choice(list(pca_bin[current_bin[0]][current_bin[1]]))
            while perm_ind not in pheno_dict or pheno_dict[perm_ind] is None:
                perm_ind = random.choice(list(pca_bin[current_bin[0]][current_bin[1]]))
            #need to make sure we match the same number as in true_phenos, might need to this several times
            #while perm_ind_pair[0] not in sample_dict or perm_ind_pair[1] not in sample_dict:
            #perm_pheno.append(pheno_dict[sample_dict[perm_ind_pair[0]]])
            #perm_pheno.append(pheno_dict[sample_dict[perm_ind_pair[1]]])
            perm_pheno.append(pheno_dict[perm_ind])
        mean_perms.append(np.mean(perm_pheno))
    return mean_perms


print 'Running permutations! [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
permutations=100
out = open(args.out, 'w')
out.write('\t'.join(['pos', 'fdr', 'num_pairs_inds']) + '\n')
for my_pos in pos_array:
    #true_phenos = pos_ibd_phenos[my_pos] #list of true phenotypes, make sure none and None
    #true_phenos_none = filter(None, true_phenos) #note: if there are zeros they will be removed
    #true_inds = list(pos_ibd_inds[my_pos]) #list of true ind tuples (randomly ordered list but tuples are right)
    true_inds = reduce(lambda x, y: x.union(y), pos_ibd_inds[my_pos], set())
    #try: this will be complicated because it's grabbing in a string formatted loop. might need to split this out
    true_phenos = []
    for my_ind in list(true_inds):
        true_phenos.append(pheno_dict[my_ind])
    
    true_bins = []
    #for my_inds in (sample_dict[true_inds[0]], sample_dict[true_inds[1]]): #list of bins (randomly ordered)
    
    
    for my_inds in true_inds: #list of bins (randomly ordered)
        true_bins.append(ind_pca_bin[my_inds])
    
    #if my_pos == 16855618:
    #    print true_phenos
    #    print true_bins
    #    print len(true_phenos)
    #    # print len(true_bins)
    #else:
    #    print len(true_phenos) #why are true_phenos and true_inds different lengths?
    #    #print len(true_inds)
    #    print true_phenos[0]
    #    #print true_inds[0]
    #    print my_pos
    #    
    #print len(true_bins)
    #print len(true_inds)
    #print true_inds[0]
    #    
    """
    in actuality (after mpg), compare true phenos to permuted with ranksum test, then ask how different z's are from 0
    i think this will provide magnitude of haplotype effect size, rather than just empirical p
    currently multiplying by 2 since there are 2 ways a distribution can outly the permuted
    """
    mean_perms = []
    mean_perms = run_perms_pca(permutations, mean_perms)
    #print mean_perms
    #if None in mean_perms:
    #    print 'mean_perms has None'
    #elif None in true_phenos:
    #    print 'true_phenos has None'
    fdr = 2*min(float(sum(np.mean(true_phenos)>mean_perms))/(len(mean_perms)+1),
              float(sum(np.mean(true_phenos)<mean_perms))/(len(mean_perms)+1))
    if fdr < 0.05:
        print 'strong assoc'
        mean_perms = run_perms_pca(900, mean_perms)
        fdr = 2*min(float(sum(np.mean(true_phenos)>mean_perms))/(len(mean_perms)+1),
              float(sum(np.mean(true_phenos)<mean_perms))/(len(mean_perms)+1))
        print [my_pos, fdr, len(true_inds)]
    out.write('\t'.join(map(str, [my_pos, fdr, len(true_inds)])) + '\n')
    out.flush()
print 'Done! [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'

out.close()
