import argparse
import collections
from geopy.distance import vincenty
import gzip
from datetime import datetime

#calculates the distance between two individuals
def calculate_distance(ind1, ind2, birth_dict, lat_lon):
    pass

def main(args):
    """
    birth_header - header for birth info (municipality, etc)
    birth_dict - dictionary containing info for every individual
    lat_lon_dict - dictionary containing codes for municipality to lat/lon
    """
    print 'starting main [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'
    match = gzip.open(args.match)
    birth_info = open(args.birth_info)
    lat_lon = open(args.lat_lon)
    header = birth_info.readline().strip().split(',')
    birth_header = {}
    for i in range(len(header)):
        birth_header[header[i]] = i
    birth_dict = {}
    for line in birth_info:
        line = line.strip().split(',')
        birth_dict[line[birth_header['ID2']]] = line
    print 'read birth record data [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'
    lat_lon_dict = {}
    lat_lon_header = {}
    header = lat_lon.readline().strip().split()
    for i in range(len(header)):
        lat_lon_header[header[i]] = i
    for line in lat_lon:
        line = line.strip().split()
        lat_lon_dict[line[lat_lon_header['CODE']]] = line
    print 'read lat/lon data [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'
    
    cum_ibd = {}
    pair_dist = {}
    i=0
    for line in match:
        i+=1
        line = line.strip().split()
        id1 = line[0]
        id2 = line[2]
        ind_pairs = tuple(sorted((id1, id2)))
        if ind_pairs in cum_ibd:
            cum_ibd[ind_pairs] += float(line[10])
        else:
            cum_ibd[ind_pairs] = float(line[10])
        if id1 in birth_dict and id2 in birth_dict:# and birth_dict[id1]['Birth.records.avail'] == '1' and birth_dict[id2]['Birth.records.avail'] == '1' and ind_pairs not in pair_dist:
            mun1 = lat_lon_dict[birth_dict[id1][birth_header['SKUNTA']]]
            ind1 = (float(mun1[lat_lon_header['LAT']]), float(mun1[lat_lon_header['LON']]))
            mun2 = lat_lon_dict[birth_dict[id2][birth_header['SKUNTA']]]
            ind2 = (float(mun2[lat_lon_header['LAT']]), float(mun2[lat_lon_header['LON']]))
            dist = vincenty(ind1, ind2).kilometers
            pair_dist[ind_pairs] = dist
        if i%10000000 == 0:
            print 'line ' + str(i) + ' [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'
            print 'cum_ibd keys: ' + str(len(cum_ibd.keys()))
            print 'pair_dist keys: ' + str(len(pair_dist.keys()))
    
    print 'start writing [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']'
    out = gzip.open(args.out, 'w')
    for inds in cum_ibd.keys():
        try:
            out.write('\t'.join(inds) + '\t' + pair_dist[inds] + '\n')
        except KeyError:
            pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--match')
    parser.add_argument('--birth_info')
    parser.add_argument('--lat_lon')
    parser.add_argument('--out')
    
    args = parser.parse_args()
    main(args)
