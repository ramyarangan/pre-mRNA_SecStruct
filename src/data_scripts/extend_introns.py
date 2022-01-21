## Last modified 1/9/2022
## Example usage: 
##  python data_scripts/extend_introns.py ../database/introns/standard_allsize_min_50_max_600/base_info.dat ../database/introns/standard_allsize_min_50_max_600_extend20/ 20

import sys
import os
from config import DATABASE_PATH
from util.gene_file_io import *

input_file = sys.argv[1] # Must be a .dat file.
output_dir = sys.argv[2]
window = int(sys.argv[3])

# Make sure the bed file has unique entries 
def get_unique_bed_data(bed_data, bps):
	new_bed_data = []
	new_bps = []
	chrpos_set = set()
	for ii, bed_data_item in enumerate(bed_data):
		if tuple(bed_data_item[:3]) not in chrpos_set:
			new_bed_data += [bed_data_item]
			new_bps += [bps[ii]]
			chrpos_set.add(bed_data_item[:3])
	return new_bed_data, new_bps

# Convert the .dat file to a .bed file if needed
base_data = read_base_data(input_file)
bed_data = [x[0][1] for x in base_data]
bps = [x[0][0] for x in base_data]

# Various duplicates in the "decoy" list need to be removed
bed_data, bps = get_unique_bed_data(bed_data, bps) 
input_bedfile = input_file.replace('.dat', '.bed')
write_bed(bed_data, filename=input_bedfile)

# Create a new bedfile with shifted genome coordinates
f = open(input_bedfile)
input_lines = f.readlines()
f.close()

output_bedfile = os.path.join(output_dir, 'base_info.bed')
output_datfile = os.path.join(output_dir, 'base_info.dat')

f = open(output_bedfile, 'w')

for input_line in input_lines:
	input_items = input_line.split('\t')
	input_items[1] = str(int(input_items[1]) - window)
	input_items[2] = str(int(input_items[2]) + window)
	input_items[3] = input_items[3].replace('\n', '')
	f.write('%s' % '\t'.join(input_items))

f.close()

genome_file = DATABASE_PATH + 'genome/sacCer3.fa'
fasta_file = fasta_from_bed(output_bedfile, genome_file)
fasta_seq_tags = read_fasta(fasta_file)
fasta_seqs = [x[1] for x in fasta_seq_tags]
clear_tmp_files()
if len(fasta_seqs) != len(bed_data):
	raise RuntimeError("fasta sequences retrieved from genome not equal to the length of the bed file")

new_base_data = []
for ii, bp in enumerate(bps):
	new_bp = str(int(bp) + window)
	new_bed_data = list(bed_data[ii]) + [str(window), str(window)]
	new_bed_data[1] = str(int(new_bed_data[1]) - window)
	new_bed_data[2] = str(int(new_bed_data[2]) + window)
	new_base_data += [((new_bp, tuple(new_bed_data)), fasta_seqs[ii])]
write_base_data(output_datfile, new_base_data)

