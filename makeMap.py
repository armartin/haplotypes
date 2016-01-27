__author__ = 'armartin'
import argparse
import sys

def main(args):
    genmap = open(args.genmap)
    bim = open(args.bim)
    chr = args.chr
    map_bim = args.map_bim
    if args.out is not None:
        my_map = open(args.out, 'w')
    elif map_bim != 'bim':
        my_map = open(args.bim.replace('bim', 'map'), 'w')
    else:
        my_map = open(args.bim + '2', 'w')
    for snp in full_map(chr, genmap, bim, map_bim):
        if map_bim != 'map':
            to_write = [snp['chr'], snp['rsid'], snp['gen_pos'], snp['phys_pos'], snp['a0'], snp['a1']]
        else:
            to_write = [snp['chr'], snp['rsid'], snp['gen_pos'], snp['phys_pos']]
        my_map.write('\t'.join(map(str, to_write)) + '\n')
    my_map.close()
    bim.close()
    genmap.close()

def full_map(chr, genmap, bim, map_bim=None, haps=None): #assert that map_bim and haps can't both have non-none values
    """
    Iterator that returns list including interpolated genetic positions
    [chr, rsid, float(interpolate), int(phys_pos), a0, a1]
    to do: set map_bim and write map output as well
    """
    genmap.readline()
    #position COMBINED_rate(cM/Mb) Genetic_Map(cM)
    #72765 0.1245577896 0
    start = genmap.readline().strip().split()
    (start_bp, start_cM) = (int(start[0]), float(start[2]))
    end = genmap.readline().strip().split()
    (end_bp, end_cM) = (int(end[0]), float(end[2]))
    
    print chr
    
    current_args = [None, start_bp, end_bp, None, bim, start_cM, end_cM, genmap, chr, None, None]
    
    for bim_line in bim:
        bim_line = bim_line.strip().split()
        if haps is not None:
            (rsid, phys_pos, a0, a1) = (bim_line[1], int(bim_line[2]), bim_line[3], bim_line[4]) #fix if map_bim=='map'
            other = bim_line[5:len(bim_line)]
        else:
            (rsid, phys_pos, a0, a1) = (bim_line[1], int(bim_line[3]), bim_line[4], bim_line[5]) #fix if map_bim=='map'
        while phys_pos < start_bp:
            proportion = (float(phys_pos) * float(start_cM)) / float(start_bp)
            if proportion < 1e-4: proportion = 1e-4
            bim_line = bim.readline().strip().split()
            if haps is not None:
                #yield final_checks([chr, rsid, str(proportion), phys_pos, a0, a1, other]) ##ISSUE IS HERE
                yield {'chr': chr, 'rsid': rsid, 'gen_pos': str(proportion), 'phys_pos': phys_pos, 'a0': a0, 'a1': a1, 'haps': other} ##ISSUE IS HERE
                (rsid, phys_pos, a0, a1) = (bim_line[1], int(bim_line[2]), bim_line[3], bim_line[4]) #fix if map_bim=='map'
            else:
                #yield final_checks([chr, rsid, str(proportion), phys_pos, a0, a1])
                yield {'chr': chr, 'rsid': rsid, 'gen_pos': str(proportion), 'phys_pos': phys_pos, 'a0': a0, 'a1': a1}
                #yield [chr, rsid, str(proportion), phys_pos, a0, a1]
                (rsid, phys_pos, a0, a1) = (bim_line[1], int(bim_line[3]), bim_line[4], bim_line[5]) #fix if map_bim=='map'
            
        (current_args[3], current_args[0], current_args[-2], current_args[-1]) = (rsid, phys_pos, a0, a1)
        while True:
            #check all conditions, break if reached the end of genetic map or bim file changed but genetic map file didn't
            new_args, to_write = check_conditions(current_args)
            
            if to_write is not None:
                if to_write['gen_pos'] < 1e-4:
                    to_write['gen_pos'] = 1e-4
                if haps is not None:
                    to_write['haps'] = other
                yield to_write
                #yield final_checks(to_write)
                break
            if new_args is None:
                break
            else:
                current_args = new_args

#don't start genetic positions quite at 0 because this throws off program (e.g. hapi-ur) assumptions
## return genetic positions for every element
def final_checks(write_vars):
    """
    fix 0 genetic positions so hapi-ur doesn't crash
    """
    print write_vars
    try:
        if float(write_vars['gen_pos']) < 1e-4: #note: other chromosomes (not 8, 11, 14, 19) don't seem to get past this point... ?
            write_vars['gen_pos'] = 1e-4
    except TypeError:
        if float(write_vars[2]) < 1e-4: #note: other chromosomes (not 8, 11, 14, 19) don't seem to get past this point... ?
            write_vars[2] = 1e-4
    return(write_vars)

def check_conditions(all_args):
    """
    this function checks where the bim file position is with respect to the genetic map
    returns same arguments that will be present next iteration and potentially what it will write - None otherwise
    """
    phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1 = all_args
    to_write = {'chr': chr,
                    'rsid': rsid,
                    'phys_pos': phys_pos,
                    'a0': a0,
                    'a1': a1}
    if phys_pos > end_bp:
        #Criteria 1 - genotypes ahead of genetic map
        while phys_pos > end_bp:
            (start_bp, start_cM) = (end_bp, end_cM)
            end = genmap.readline().strip().split() 
            if end == []:
                proportion = (phys_pos * end_cM) / end_bp
                #to_write = [chr, rsid, str(proportion), phys_pos, a0, a1]
                to_write['gen_pos'] = phys_pos
                return ([phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1], to_write)
            else:
                (end_bp, end_cM) = (int(end[0]), float(end[2]))
        return ([phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1], None)
    else:
        #Criteria 2 - genotypes not ahead of genetic map
        if phys_pos == start_bp:
            #to_write = [chr, rsid, start_cM, phys_pos, a0, a1]
            to_write['gen_pos'] = start_cM
            return ([phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1], to_write)
        elif phys_pos == end_bp:
            #to_write = [chr, rsid, end_cM, phys_pos, a0, a1]
            to_write['gen_pos'] = end_cM
            return ([phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, chr, a0, a1], to_write)
        elif phys_pos > start_bp and phys_pos < end_bp:
            interpolate = (phys_pos - start_bp)*(end_cM - start_cM)/(end_bp - start_bp) + start_cM
            #to_write = [chr, rsid, interpolate, phys_pos, a0, a1]
            to_write['gen_pos'] = interpolate
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