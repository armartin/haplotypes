import argparse
import gzip

def open_file(filename):
    if filename.endswith('gz'):
        my_file = gzip.open(filename)
    else:
        my_file = open(filename)
    return my_file
    
def main(args):
    gaps = open_file(args.gaps)
    match = open_file(args.match)
    out = gzip.open(args.out, 'w')
    
    gaps.readline()
    gap_starts = {}
    gap_ends = {}
    for line in gaps:
        line = line.strip().split()
        if line[1] in gap_starts:
            gap_starts[line[1]].append(int(line[2]))
            gap_ends[line[1]].append(int(line[3]))
        else:
            gap_starts[line[1]] = [int(line[2])]
            gap_ends[line[1]] = [int(line[3])]
    
    for line in match:
        line = line.strip().split()
        chr = int(line[4])
        start = int(line[5])
        end = int(line[6])
        in_gap = False
        for gap in range(len(gap_starts[chr])):
            if start > gap_starts[chr][gap] and start < gap_ends[chr][gap] or end > gap_starts[chr][gap] and end < gap_ends[chr][gap] or start < gap_starts[chr][gap] and end > gap_starts[chr][gap]:
                in_gap = True
                continue
        if not in_gap:
            out.write('\t'.join(line) + '\n')
    out.close()

if __name__ == '__main__':    
    parser = argparse.ArgumentParser()
    parser.add_argument('--gaps', default='/home/unix/armartin/finrisk/haplotypes/b37_gaps.txt')
    parser.add_argument('--match', required=True)
    parser.add_argument('--out', required=True)
    
    args = parser.parse_args()
    main(args)