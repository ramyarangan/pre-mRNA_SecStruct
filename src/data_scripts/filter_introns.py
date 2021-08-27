import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]

input_file = open(input_filename)
input_lines = input_file.readlines()
input_file.close()

LEN_MAX = 600
LEN_MIN = 50

output_lines = []
for ii in range(0, int(len(input_lines)/2)):
	bp_pos_line = input_lines[ii * 2]
	bp_pos = bp_pos_line.split()[0]
	intron = input_lines[ii * 2 + 1]

	if (int(bp_pos) == -1):
		continue

	if (len(intron) > LEN_MAX) or \
		(len(intron) < LEN_MIN):
		continue

	output_lines = output_lines + [bp_pos_line, intron]

output_file = open(output_filename, 'w')

for output_line in output_lines:
	output_file.write(output_line)

output_file.close()