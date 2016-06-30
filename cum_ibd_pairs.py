import argparse
import gzip
import collections
from datetime import datetime
from geopy.distance import vincenty


## gets the latitude/longitude for an individual
def ind_loc(lat_lon_dict, lat_lon_header, birth_dict, birth_header, ind_id):
    try:
        mun = lat_lon_dict[birth_dict[ind_id][birth_header['SKUNTA']]]
        #reg = lat_lon_dict[birth_dict[ind_id][birth_header['SKUNTA_LAAN']]]
        ind = (float(mun[lat_lon_header['LAT']]), float(mun[lat_lon_header['LON']]))
        return(ind)
    except KeyError:
        return 'NA'

def main(args):
    if args.match.endswith('gz'):
        match = gzip.open(args.match)
    else:
        match = open(args.match)
    
    ## first, make sure all pairs are represented
    print 'Getting all pairs [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'
    inds = open(args.inds)
    #inds.readline()
    #inds.readline()
    samples = []
    for line in inds:
        line = line.strip().split()
        samples.append(line[1])
    
    cum_ibd = {}
    all_ibd = {}
    for ind1 in samples:
        for ind2 in samples:
            if ind1 != ind2:
                ind_pairs = tuple(sorted((ind1, ind2)))
                cum_ibd[ind_pairs] = 0
                all_ibd[ind_pairs] = []
    
    ## add birth locs where possible
    if args.birth is not None and args.lat_lon is not None:
        print 'Adding birth locations [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'
        birth = open(args.birth)
        header = birth.readline().strip().split(',')
        birth_header = {}
        for i in range(len(header)):
            birth_header[header[i]] = i
        birth_dict = {}
        for line in birth:
            line = line.strip().split(',')
            birth_dict[line[birth_header['ID2']]] = line
        
        lat_lon = open(args.lat_lon)
        lat_lon_dict = {}
        lat_lon_header = {}
        header = lat_lon.readline().strip().split()
        for i in range(len(header)):
            lat_lon_header[header[i]] = i
        for line in lat_lon:
            line = line.strip().split()
            lat_lon_dict[line[lat_lon_header['CODE']]] = line
    
    ## next, count cumulative IBD sharing across pairs and distance between them
    i=0
    pair_dist = {}
    for line in match:
        if i%10000000 == 0:
            print 'line ' + str(i) + ' [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'
        i+=1
        line = line.strip().split()
        ind1 = line[1].split('.')[0]
        ind2 = line[3].split('.')[0]
        if ind1 != ind2:
            ind_pairs = tuple(sorted((ind1, ind2)))
            cum_ibd[ind_pairs] += float(line[10])
            all_ibd[ind_pairs].append(float(line[10]))
            if args.birth is not None and args.lat_lon is not None:
                ind1_loc = ind_loc(lat_lon_dict, lat_lon_header, birth_dict, birth_header, ind1)
                ind2_loc = ind_loc(lat_lon_dict, lat_lon_header, birth_dict, birth_header, ind2)
                if ind1_loc != 'NA' and ind2_loc !='NA':
                    dist = vincenty(ind1_loc, ind2_loc).kilometers
                    pair_dist[ind_pairs] = dist
                    if i%10000000 == 0:
                        print ind_pairs
                        print pair_dist[ind_pairs]
    
    ## write cumulative IBD
    print 'Writing output [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'
    all_pairs = cum_ibd.keys()
    out = gzip.open(args.out, 'w')
    for pair in all_pairs:
        if args.birth is not None:
            if pair in pair_dist:
                out.write(pair[0] + '\t' + pair[1] + '\t' + str(cum_ibd[pair]) + '\t' + ','.join(map(str, sorted(all_ibd[pair]))) + '\t' + str(pair_dist[pair]) + '\t' +
                          birth_dict[pair[0]][birth_header['SKUNTA_LAAN']] + '\t' + birth_dict[pair[1]][birth_header['SKUNTA_LAAN']] + '\n')
            else:
                out.write(pair[0] + '\t' + pair[1] + '\t' + str(cum_ibd[pair]) + '\t' + ','.join(map(str, sorted(all_ibd[pair]))) + '\tNA\n')
        else:
            out.write(pair[0] + '\t' + pair[1] + '\t' + str(cum_ibd[pair]) + '\t' + ','.join(map(str, sorted(all_ibd[pair]))) + '\n')
    out.close()


if __name__ == '__main__':    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--match')
    parser.add_argument('--inds')
    parser.add_argument('--birth')
    parser.add_argument('--lat_lon')
    parser.add_argument('--out')
    
    args = parser.parse_args()
    main(args)