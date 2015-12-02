__author__ = 'armartin'
import argparse
import sys

"""
parser.add_argument('--chr', default='22') #will replace chr[\d] with chr[chr]
parser.add_argument('--genmap')
parser.add_argument('--bim')
parser.add_argument('--map_bim', default='bim') #read/write map/bim files
parser.add_argument('--out')
"""
def main(args):
    genmap = open(args.genmap)
    bim = open(args.bim)
    chr = args.chr
    map_bim = args.map_bim
    if args.out is not None:
        my_map = open(args.out, 'w')
    else:
        my_map = open(args.bim.replace('bim', 'map'), 'w')
    
    for snp in full_map(chr, genmap, bim, map_bim):
        my_map.write('\t'.join(map(str, snp)) + '\n')
    my_map.close()
    bim.close()
    genmap.close()

def full_map(chr, genmap, bim, map_bim):
    """
    reads through first line of genetic and 
    """
    genmap.readline()
    #position COMBINED_rate(cM/Mb) Genetic_Map(cM)
    #72765 0.1245577896 0
    start = genmap.readline().strip().split()
    (start_bp, start_cM) = (start[0], start[2])
    end = genmap.readline().strip().split()
    (end_bp, end_cM) = (end[0], end[2])
    
    print chr
    
    
    
    for bim_line in bim:
        bim_line = bim_line.strip().split()
        (rsid, phys_pos, a0, a1) = (bim_line[1], bim_line[3], bim_line[4], bim_line[5])
        while phys_pos < start_bp:
            proportion = (float(phys_pos) * float(start_cM)) / float(start_bp)
            yield final_checks([chr, rsid, str(proportion), phys_pos, a0, a1])
            bim_line = bim.readline().strip().split()
            (rsid, phys_pos, a0, a1) = (bim_line[1], bim_line[3], bim_line[4], bim_line[5])
        print 'paragraph 1'
        current_args = [phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1]
        while True:
            current_args, to_write = check_conditions(current_args)
            if to_write is not None:
                yield final_checks(to_write)
            if current_args is None:
                break

#don't start genetic positions quite at 0 because this throws off program (e.g. hapi-ur) assumptions
## return genetic positions for every element
def final_checks(write_vars):
    write_vars[2] = str(write_vars[2])
    if str(write_vars[2]) == '0.0':
        write_vars[2] = 1e-4
    return(write_vars)

"""
this function checks where the bim file position is with respect to the genetic map
returns same arguments that will be present next iteration and potentially what it will write - None otherwise
"""
def check_conditions(all_args):
    phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1 = all_args
    if int(phys_pos) > int(end_bp):
        #Criteria 1 - genotypes ahead of genetic map
        while int(phys_pos) > int(end_bp):
            (start_bp, start_cM) = (end_bp, end_cM)
            end = genmap.readline().strip().split() 
            if end == []:
                proportion = (float(phys_pos) * float(end_cM)) / float(end_bp)
                to_write = [chr, rsid, str(proportion), phys_pos, a0, a1]
                return ([phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1], to_write)
            else:
                (end_bp, end_cM) = (end[0], end[2])
        return ([phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1], None)
    else:
        #Criteria 2 - genotypes not ahead of genetic map
        if phys_pos == start_bp:
            to_write = [chr, rsid, str(start_cM), phys_pos, a0, a1]
            return ([phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1], to_write)
        elif phys_pos == end_bp:
            to_write = [chr, rsid, str(end_cM), phys_pos, a0, a1]
            return ([phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1], to_write)
        elif int(phys_pos) > int(start_bp) and int(phys_pos) < int(end_bp):
            interpolate = float((int(phys_pos) - int(start_bp))*(float(end_cM) - float(start_cM)))/(int(end_bp) - int(start_bp)) + float(start_cM)
            to_write = [chr, rsid, str(interpolate), phys_pos, a0, a1]
            return ([phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1], to_write)
        else:
            print 'Criteria 2 - Something wrong when genetic data is before physical positions'
        return (None, None)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse some args')
    parser.add_argument('--chr', default='22') #will replace chr[\d] with chr[chr]
    parser.add_argument('--genmap')
    parser.add_argument('--bim')
    parser.add_argument('--map_bim', default='bim') #read/write map/bim files
    parser.add_argument('--out')
    
    args = parser.parse_args()
    main(args)