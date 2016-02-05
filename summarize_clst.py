import argparse
import gzip

## open gzipped and plain files
def open_file(filename):
    if filename.endswith('gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return my_file

def main(args):
    dash = open_file(args.dash)
    out = gzip.open(args.out, 'w')
    if args.pheno is not None:
        pheno = open_file(args.pheno)
        pheno_inds = set()
        for line in pheno:
            line = line.strip().split()
            pheno_inds.add(line[1])
        
    for line in dash:
        line = line.strip().split()
        if args.pheno is not None:
            pheno_num = 0
            cluster_inds = line[5:len(line):2]
            for ind in cluster_inds:
                if ind in pheno_inds:
                    pheno_num += 1
            out.write('\t'.join([args.chr + '_' + line[0], line[1], line[2], str(pheno_num)]) + '\n')
        else:
            out.write('\t'.join([args.chr + '_' + line[0], line[1], line[2], str((len(line) - 5)/2)]) + '\n')
    out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dash', required=True)
    parser.add_argument('--chr', required=True)
    parser.add_argument('--pheno')
    parser.add_argument('--out')#, required=True)
    
    args = parser.parse_args()
    main(args)