import argparse
import collections
from geopy.distance import vincenty
import gzip

#calculates the distance between two individuals
def calculate_distance(ind1, ind2, birth_dict, lat_lon):
    pass

def main(args):
    """
    birth_header - header for birth info (municipality, etc)
    birth_dict - dictionary containing info for every individual
    lat_lon_dict - dictionary containing codes for municipality to lat/lon
    """
    match = gzip.open(args.match)
    birth_info = open(args.birth_info)
    lat_lon = open(args.lat_lon)
    header = birth_info.readline().strip().split(',')
    birth_header = {}
    for i in range(len(header)):
        birth_header[header[i]] = i
    birth_dict = {}
    for line in birth_info:
        line = line.strip().split()
        birth_dict[birth_header['ID2']] = line
    lat_lon_dict = {}
    lat_lon_header = {}
    header = lat_lon.readline().strip().split()
    for i in range(len(header)):
        lat_lon_header[header[i]] = i
    for line in lat_lon:
        line = line.strip().split()
        lat_lon_dict[lat_lon_header['CODE']] = line
    
    cum_ibd = {}
    pair_dist = {}
    for line in match:
        line = line.strip().split()
        id1 = line[0]
        id2 = line[2]
        ind_pairs = sorted([id1, id2])
        print ind_pairs
        if ind_pairs in cum_ibd:
            cum_ibd[ind_pairs] += float(line[11])
        else:
            cum_ibd[ind_pairs] = float(line[11])
        if id1 in birth_dict and id2 in birth_dict and birth_dict[id1]['Birth.records.avail'] == '1' and birth_dict[id2]['Birth.records.avail'] == '1':
            mun1 = lat_lon_dict[birth_dict[id1]['SKUNTA']]
            ind1 = (float(mun1[lat_lon_header['LAT']]), float(mun1[lat_lon_header['LON']]))
            mun2 = lat_lon_dict[birth_dict[id2]['SKUNTA']]
            ind2 = (float(mun2[lat_lon_header['LAT']]), float(mun2[lat_lon_header['LON']]))
            dist = vincenty(ind1, ind2).kilometers
            pair_dist[ind_pairs] = dist
    
    out = gzip.open(args.out, 'w')
    for inds in cum_ibd.keys():
        out.write('\t'.join(inds) + '\t' + pair_dist[inds] + '\n')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--match')
    parser.add_argument('--birth_info')
    parser.add_argument('--lat_lon')
    parser.add_argument('--out')
    
    args = parser.parse_args()
    main(args)
