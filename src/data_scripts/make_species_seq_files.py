"""
Make evolutionarily related intron sets

Includes functionality to make intron control sets that match the 
mutation / in-del level compared to S. cerevisiae.

Example function call:

python data_scripts/make_species_seq_files.py ../database/hooks_alignments/ species_hooks --min 50 --max 600 --make_shuffle --make_shuffle_seq_matched --make_phylo
"""


import argparse
import os
import random 
import numpy as np

from config import DATABASE_PATH
from make_case_matched_controls import make_shuffled_control
from make_case_matched_controls import make_shuffle_seq_matched_control
from util.gene_file_io import * 
from util.aln_util import *

def get_seq(seq, pos_info):
	[five_ss, bp, three_ss] = pos_info
	full_seq = seq[five_ss:(three_ss+1)].upper()
	return full_seq

def make_mutations(scer_intron, scer_intron_exp, species_seq):
	# Tabulate modification counts
	insertions = 0
	deletions = 0
	mutations = 0
	for ii in range(len(scer_intron_exp)):
		if (scer_intron_exp[ii] == ".") or (scer_intron_exp[ii] == "-"):
			if (species_seq[ii] != ".") and (species_seq[ii] != "-"):
				insertions += 1
		elif (species_seq[ii] == ".") or (species_seq[ii] == "-"):
			if (scer_intron_exp[ii] != ".") or (scer_intron_exp[ii] != "-"):
				deletions += 1
		elif (scer_intron_exp[ii] != species_seq[ii]):
			mutations += 1
	print(str(insertions) + " " + str(deletions) + " " + str(mutations))
	# Make a np array of characters
	intron = np.array(list(scer_intron))
	
	# Do the deletions
	keep_locs = random.sample(range(len(scer_intron)), len(scer_intron) - deletions)
	keep_locs.sort()
	intron = intron[keep_locs]
	
	# Make space for the insertions
	insert_intron = np.array([" "] * (len(intron) + insertions))
	intron_locs = random.sample(range(len(insert_intron)), len(intron))
	intron_locs.sort()
	insert_intron[intron_locs] = intron
	intron = insert_intron

	# Do the mutations and insertions
	for ii in range(len(intron)):
		if intron[ii] == " ":
			intron[ii] = get_rnd_char()

	mut_locs = random.sample(range(len(intron)), mutations)
	for mut_loc in mut_locs:
		old_char = intron[mut_loc]
		while old_char == intron[mut_loc]:
			intron[mut_loc] = get_rnd_char()

	# Return in string form
	return "".join(intron)



# First make species: [(Gene name, intron seq, branchpoint pos, intron num), ...] dictionary.
def make_control_seqs(alignment_dir):
	control_seqs = {}

	for filename in os.listdir(alignment_dir):
		if filename.endswith(".stk"):
			f = open(alignment_dir + filename)
			lines = f.readlines()
			f.close()

			gene_name = filename.split(".")[0]

			annotation_str = lines[len(lines)-2].split()[2]
			intron_pos_list = get_introns(annotation_str)
			print(filename)
			(num_skipped, starting_pos) = get_skips(lines)

			for ii in range(len(intron_pos_list)):
				intron_pos = intron_pos_list[ii]

				scer_intron_exp = ""
				scer_intron = "" 
				for jj in range(len(lines)-num_skipped):
					line = lines[jj + starting_pos]
					species_name = (line.split()[0]).split('_')[0]
					if (species_name == "scer"):
						scer_intron_exp = get_seq(line.split()[1], intron_pos)
						[scer_intron, _] = condense_seq(line.split()[1], intron_pos)

				for jj in range(len(lines)-num_skipped):
					line = lines[jj + starting_pos]
					species_name = (line.split()[0]).split('_')[0]
					species_seq = get_seq(line.split()[1], intron_pos)
					[condensed_seq, bp] = condense_seq(line.split()[1], intron_pos)
					print(species_name)
					if len(condensed_seq) == 0:
						continue
					mutated_seq = make_mutations(scer_intron, scer_intron_exp, species_seq)
					# print(len(mutated_seq)) # Useful check that mutations are reasonable
					# print(len(condensed_seq))

					new_entry = (gene_name, mutated_seq, bp, ii + 1)
					if species_name in control_seqs:
						control_seqs[species_name] += [new_entry]
					else:
						control_seqs[species_name] = [new_entry]

	return control_seqs


def write_seqdat_files(species_seqs, seq_dir, min_size=-1, max_size=-1):
	# Get number of introns for S. cer for labeling. 
	intron_cnts = {}
	for intron in species_seqs['scer']:
		gene_name = intron[0]
		if gene_name in intron_cnts:
			intron_cnts[gene_name] = max(intron_cnts[gene_name], intron[3])
		else:
			intron_cnts[gene_name] = 1

	for species_name, introns in species_seqs.items():
		if not os.path.exists(seq_dir):
			os.mkdir(seq_dir)
		dir_path = os.path.join(seq_dir, species_name + '/')
		if not os.path.exists(dir_path):
			os.mkdir(dir_path)
		filename = os.path.join(dir_path, 'base_info.dat')
		
		f = open(filename, 'w')

		for intron in introns:
			gene_name = intron[0]
			if gene_name not in intron_cnts:
				continue

			if intron_cnts[gene_name] > 1:
				gene_name += '_' + str(intron[3])

			seq = intron[1]

			# Branchpoint -1
			if int(intron[2]) == -1:
				continue

			if (min_size != -1) and \
				((len(seq) < min_size) or (len(seq) > max_size)):
				continue
			end_pos = 1 + len(seq)
			# Make a fake chromosome location for now rather than finding
			# real genome position in genome annotations
			info_str = str(intron[2]) + ' chrI\t1\t' + str(end_pos) + \
				'\t\t0\t+'
			f.write('%s\n' % info_str)
			f.write('%s\n' % intron[1])

		f.close()

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Parameters for processing intron alignment data into intron classes')
	parser.add_argument('alignment_dir', type=str, help='Path to directory storing alignments in Stockholm format')
	parser.add_argument('species_intron_class_dir', type=str, help='Path that will store base_info.dat files for all species')
	parser.add_argument('--make_shuffle', default=False, action='store_true', \
		 help='Make a shuffled intron case-matched control set')
	parser.add_argument('--make_shuffle_seq_matched', default=False, action='store_true', \
		 help='Make a shuffle intron case-matched control set with matching splicing sequences')
	parser.add_argument('--make_phylo', default=False, action='store_true', \
		 help='Make a control set with mutation and indel frequencies from Scer matched')
	parser.add_argument('--min', default=-1, type=int, help='Minimum intron length')
	parser.add_argument('--max', default=-1, type=int, help='Maximum intron length')
	args = parser.parse_args()

	alignment_dir = args.alignment_dir
	species_intron_class_dir = args.species_intron_class_dir
	make_shuffle = args.make_shuffle
	make_shuffle_seq_matched = args.make_shuffle_seq_matched
	make_phylo = args.make_phylo
	min_size = args.min
	max_size = args.max

	species_seqs = make_species_seqs(alignment_dir)
	file_ext = ""
	if min_size != -1:
		file_ext = "_min_" + str(min_size) + "_max_" + str(max_size)

	standard_dir = os.path.join(DATABASE_PATH, \
		'introns/' + species_intron_class_dir + '/standard' + file_ext + '/')
	write_seqdat_files(species_seqs, standard_dir, min_size=min_size, max_size=max_size)

	if make_phylo:
		phylo_dir = species_intron_class_dir + '/standard' + file_ext + '_phylo_control/'
		phylo_dir = os.path.join(DATABASE_PATH, 'introns/' + phylo_dir)
		control_seqs = make_control_seqs(alignment_dir)
		write_seqdat_files(control_seqs, phylo_dir, min_size=min_size, max_size=max_size)

	if make_shuffle:
		for species_name, _ in species_seqs.items():
			dir_path = os.path.join(standard_dir, species_name + '/')
			filename = os.path.join(dir_path, 'base_info.dat')		
			base_data = read_base_data(filename)
			intron_class = species_intron_class_dir + '/standard' + file_ext
			make_shuffled_control(base_data, intron_class, \
				do_species=True, species_name=species_name)

	if make_shuffle_seq_matched:
		for species_name, _ in species_seqs.items():
			dir_path = os.path.join(standard_dir, species_name + '/')
			filename = os.path.join(dir_path, 'base_info.dat')		
			base_data = read_base_data(filename)
			intron_class = species_intron_class_dir + '/standard' + file_ext
			make_shuffle_seq_matched_control(base_data, intron_class, \
				do_species=True, species_name=species_name)

