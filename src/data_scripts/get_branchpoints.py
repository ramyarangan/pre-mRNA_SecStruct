import sys

# For now only look for consensus branchpoint sequences
consensus = 'TACTAAC'

fasta_filename = sys.argv[1]
bed_filename = sys.argv[2]
out_filename = sys.argv[3]

fasta_file = open(fasta_filename)
fasta_lines = fasta_file.readlines()
fasta_file.close()

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

no_bp_cnt = 0
bp_cnt_rev = 0
bp_cnt_for = 0
bp_lines = []
for ii in range(0, len(bed_lines)):
	bed_items = bed_lines[ii].split()
	is_pos = (bed_items[len(bed_items)-1] == '+')
	next_str = fasta_lines[2 * ii + 1].upper()
	if not is_pos:
		next_str = reverse_invert(next_str)
	if (next_str[len(next_str) - 1] == '\n'):
		next_str = next_str[:-1]
	bp_idx = next_str.find(consensus)
	if (bp_idx != -1):
		bp_idx = bp_idx + 5 # Get to the right "A" in the consensus sequence
		no_bp_cnt = no_bp_cnt + 1
	else:
		if not is_pos:
			bp_cnt_rev = bp_cnt_rev + 1
		if is_pos:
			bp_cnt_for = bp_cnt_for + 1
	bp_lines = bp_lines + [str(bp_idx) + " " + bed_lines[ii], next_str + '\n']

print(no_bp_cnt)
print(bp_cnt_rev)
print(bp_cnt_for)
print(len(bed_lines))

out_file = open(out_filename, 'w')
for bp_line in bp_lines:
	out_file.write(bp_line)
out_file.close()
