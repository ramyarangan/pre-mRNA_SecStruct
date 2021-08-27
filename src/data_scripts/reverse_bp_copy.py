import sys

# Reverse any reversed intron sequences and copy bp positions from another
# file. This is used to set up control comparisons. 
fasta_filename = sys.argv[1]
bed_filename = sys.argv[2]
bp_filename = sys.argv[3]
out_filename = sys.argv[4]

fasta_file = open(fasta_filename)
fasta_lines = fasta_file.readlines()
fasta_file.close()

other_bp_file = open(bp_filename)
other_bp_lines = other_bp_file.readlines()
other_bp_file.close()

bed_file = open(bed_filename)
bed_lines = bed_file.readlines()
bed_file.close()

def reverse_invert(intron):
	replace_ch = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	reversed_intron = ''
	for ii in range(0, len(intron)):
		if (intron[ii] == '\n'):
			continue
		reversed_intron = replace_ch[intron[ii]] + reversed_intron
	return reversed_intron

bp_lines = []
for ii in range(0, len(bed_lines)):
	bed_items = bed_lines[ii].split()
	is_pos = not (bed_items[len(bed_items)-1] == '-')
	next_str = fasta_lines[2 * ii + 1].upper()
	if not is_pos:
		next_str = reverse_invert(next_str)
	if (next_str[len(next_str) - 1] == '\n'):
		next_str = next_str[:-1]
	bp_idx = str(other_bp_lines[2 * ii])
	if (bp_idx[len(bp_idx) - 1] == '\n'):
		bp_idx = bp_idx[:-1]
	bp_lines = bp_lines + [bp_idx, next_str]

out_file = open(out_filename, 'w')
for bp_line in bp_lines:
	out_file.write(bp_line)
	out_file.write('\n')
out_file.close()
