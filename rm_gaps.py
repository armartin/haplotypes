import argparse
import gzip

def open_file(filename):
    if filename.endswith('gz'):
        gzip.open(filename)
    else:
        open(filename)
    
def main(args):
    gap = open_file(args.gap)
    match = open_file(args.match)
    out = gzip.open(args.out, 'w')
    
    gap.readline()
    gap_starts = {}
    gap_ends = {}
    for line in gap:
        line = line.strip().split()
        if line[1] in gap_starts:
            gap_starts.append(int(line[2]))
            gap_ends.append(int(line[3]))
        else:
            gap_starts = [int(line[2])]
            gap_ends = [int(line[3])]
    
    for line in match:
        line = line.strip().split()
        start = line[5]
        end = line[6]
        for gap in range(len(gap_starts)):
            if start > gap_starts[gap] and start < gap_ends[gap] or end > gap_starts[gap] and end < gap_ends[gap]:
                pass
            else:
                out.write('\t'.join(line))
    out.close()

if __name__ == '__main__':    
    parser = argparse.ArgumentParser()
    parser.add_argument('--gaps', required=True, default='/home/unix/armartin/finrisk/haplotypes/b37_gaps.txt')
    parser.add_argument('--match', required=True)
    parser.add_argument('--out', required=True)
    
    args = parser.parse_args()
    main(args)