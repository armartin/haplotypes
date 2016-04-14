import argparse
import collections
import re

parser = argparse.ArgumentParser(description='Parse some args')
parser.add_argument('--hist_files')
parser.add_argument('--out')

args = parser.parse_args()

hist_files = open(args.hist_files)
count_dict = collections.defaultdict(dict)
my_range = range(1,102)

for line in hist_files:
    line = line.strip()
    my_file = open(line)
    chr = re.search('chr(\d+)', line).group(0)
    my_file.readline()
    min_max = my_file.readline().strip()
    if 'outside' in min_max:
        outside = int(min_max.split()[1])
        my_file.readline()
    for line_count in my_file:
        line_count = line_count.strip().split()
        start = int(float(line_count[0]))
        end = int(float(line_count[2]))
        if line_count[3] == '[':
            count = int(line_count[4].replace(']:', ''))
        else:
            count = line_count[3].replace('[', '')
            count = int(count.replace(']:', ''))
        print count
        if start in count_dict:
            count_dict[start][end] = count_dict[start][end] + count
        else:
            count_dict[start][end] = count
    if 'outside' in count_dict:
        count_dict['outside'] = count_dict['outside'] + outside
    else:
        count_dict['outside'] = outside
    print count_dict

print args.out

out = open(args.out, 'w')
out.write('\t'.join(['start', 'end', 'count']) + '\n')
for i in range(len(my_range)-1):
    start = my_range[i]
    end = my_range[i+1]
    out.write('\t'.join([str(start), str(end), str(count_dict[start][end])]) + '\n')
out.write('\t'.join(['outside', 'outside', str(count_dict['outside'])]) + '\n')
out.close()
