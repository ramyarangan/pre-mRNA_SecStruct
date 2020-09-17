import sys

input_bed_filename = sys.argv[1]
output_bed_filename = sys.argv[2]
shift_dist = int(sys.argv[3])

input_bed_file = open(input_bed_filename)
bed_lines = input_bed_file.readlines()
input_bed_file.close()

control_lines = []
for bed_line in bed_lines:
	bed_items = bed_line.split()
	control_items = bed_items
	if (bed_items[len(bed_items)-1] == '-'):
		control_items[1] = str(int(control_items[1]) - shift_dist)
		control_items[2] = str(int(control_items[2]) - shift_dist)	
	else:
		control_items[1] = str(int(control_items[1]) + shift_dist)
		control_items[2] = str(int(control_items[2]) + shift_dist)
	
	control_line = '\t'.join(control_items)
	control_lines = control_lines + [control_line]

output_bed_file = open(output_bed_filename, 'w')

for line in control_lines:
	output_bed_file.write(line)
	output_bed_file.write('\n')
output_bed_file.close()