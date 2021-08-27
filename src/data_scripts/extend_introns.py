import sys

input_bedfile = sys.argv[1]
output_bedfile = sys.argv[2]
window = int(sys.argv[3])

f = open(input_bedfile)
input_lines = f.readlines()
f.close()

f = open(output_bedfile, 'w')

for input_line in input_lines:
	input_items = input_line.split('\t')
	input_items[1] = str(int(input_items[1]) - window)
	input_items[2] = str(int(input_items[2]) + window)
	input_items[3] = input_items[3].strip('\n')
	f.write('%s\n' % '\t'.join(input_items))

f.close()