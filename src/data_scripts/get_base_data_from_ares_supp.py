# 8/27/21
# Creates base_info data for an intron set from supplementary info 
# in the Talkish, et al. 2020 proto intron paper 

# Run like: 
# python data_scripts/get_base_data_from_ares_supp.py proto_allsize
# python data_scripts/get_base_data_from_ares_supp.py proto_allsize --min 50 --max 600

import pandas as pd 
import os
import argparse

from config import DATABASE_PATH
from core.gene import GeneSet
from util.gene_file_io import write_base_data

parser = argparse.ArgumentParser(description='Parameters for processing intron data')
parser.add_argument('intron_class', type=str, help='Intron class name')
parser.add_argument('--min', default=-1, type=int, help='Minimum intron length')
parser.add_argument('--max', default=-1, type=int, help='Maximum intron length')
args = parser.parse_args()

intron_class = args.intron_class
length_min = args.min
length_max = args.max

intron_data_path = 'introns/' + intron_class + "/" + intron_class + '_ares.csv'
csv_path = os.path.join(DATABASE_PATH, intron_data_path)

df = pd.read_csv(csv_path)
chr_vals = df['Chromosome'].tolist()
start_idxs = df['Start'].tolist()
end_idxs = df['End'].tolist()
strand = df['Strand'].tolist()
seqs = df['sequence'].tolist()
bp_pos = df["5'SS to BP (bp)"].tolist()

# Locate the gene for each proto intron
gene_set = GeneSet()
all_genes = gene_set.genes
genes = []
for ii, chr_val in enumerate(chr_vals):
	# Some introns are in UTR's and won't be in a gene bed coordinates
	# For those introns, just don't assign a name for now
	gene_name = ""
	for gene in all_genes: 
		if gene.chr_num != chr_val:
			continue
		(chr_start, chr_end) = gene.chr_pos
		if chr_start < start_idxs[ii] and \
			chr_end > end_idxs[ii]:
			if seqs[ii] not in gene.seq: 
				print(seqs[ii])
				continue
				# raise RuntimeError("Gene sequence does not contain intron sequence")
			gene_name = gene.refseq_name
			break
	genes += [gene_name]

base_info_items = []
base_info_seqs = []

for ii, chr_val in enumerate(chr_vals): 
	# Check that intron meets length requirements if any
	if (length_min == -1) or \
		(len(seqs[ii]) > length_min and len(seqs[ii]) < length_max):
		# sequence[bp_pos-1] gives the branchpoint A, aka TACTAAC --> position 6
		base_info_item = (bp_pos[ii] - 1, chr_vals[ii], start_idxs[ii], \
			end_idxs[ii], genes[ii], 0, strand[ii])
		base_info_items	+= [base_info_item]
		base_info_seqs += [seqs[ii]]

# Make a intron class name for this set, depending on what the length min and 
# max were that were used for filtering sequences
classname = intron_class
if length_min != -1:
	classname = intron_class + '_min_' + str(length_min) + '_max_' + str(length_max)
filepath = os.path.join(DATABASE_PATH, 'introns/' + classname + '/')
if not os.path.exists(filepath):
	os.mkdir(filepath)
base_info_path = os.path.join(filepath, 'base_info.dat')

write_base_data(base_info_path, base_info_items, base_info_seqs)

