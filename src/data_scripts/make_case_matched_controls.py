# 8/27/21
# Creates shuffled and shifted controls given an intron class
# Each intron in the base_info.dat file is line-matched to the original
# intron class (order matters in the base_info.dat files!)
# Run like: 
# python data_scripts/make_case_matched_controls.py standard_allsize_min_50_max_600 --make_shuffle --make_shifted

import os
import argparse
import random 
from config import DATABASE_PATH
from util.gene_file_io import * 

def make_shifted_control(base_data, intron_class):
	new_classname = intron_class + '_shift_' + str(shift_dist)
	new_filepath = os.path.join(DATABASE_PATH, 'introns/' + new_classname + '/')
	if not os.path.exists(new_filepath):
		os.mkdir(new_filepath)

	bed_data = []
	bps = []
	offsets = []
	for data_line in base_data: 
		bps += [data_line[0][0]]
		bed_start = int(data_line[0][1][1])
		bed_end = int(data_line[0][1][2])
		strand_dir = data_line[0][1][4]
		if strand_dir == '-':
			bed_start = bed_start - shift_dist
			bed_end = bed_end - shift_dist
		else:
			bed_start = bed_start + shift_dist
			bed_end = bed_end + shift_dist
		cur_bed_data = list(data_line[0][1])
		cur_bed_data[1] = bed_start
		cur_bed_data[2] = bed_end
		bed_data += [tuple(cur_bed_data[:5])]
		if len(cur_bed_data) > 6:
			offsets += [(cur_bed_data[5], cur_bed_data[6])]
		
	bed_filename = write_bed(bed_data)

	base_data_path = os.path.join(new_filepath, 'base_info.dat')
	write_base_data_from_bed(base_data_path, bed_filename, bps, offsets=offsets)
	clear_tmp_files()

def make_shuffled_control(base_data, intron_class, do_species=False, species_name=""):
	new_classname = intron_class + '_shuffle'
	new_filepath = os.path.join(DATABASE_PATH, 'introns/' + new_classname + '/')
	if not os.path.exists(new_filepath):
		os.mkdir(new_filepath)
	if do_species:
		new_filepath = os.path.join(new_filepath, species_name + '/')
		if not os.path.exists(new_filepath):
			os.mkdir(new_filepath)
	base_data_path = os.path.join(new_filepath, 'base_info.dat')

	base_data_items = []
	base_data_seqs = []
	for data_line in base_data: 
		((bp, bed_items), seq) = data_line
		base_data_items += [tuple([bp] + list(bed_items))]
		base_data_seqs += [''.join(random.sample(seq, len(seq)))]
	
	write_base_data_items(base_data_path, base_data_items, base_data_seqs)

def make_seq_matched(base_data, base_classname, shifted=True, \
	do_species=False, species_name=""):
	new_classname = base_classname + '_seq_matched'
	new_filepath = os.path.join(DATABASE_PATH, 'introns/' + new_classname + '/')
	if not os.path.exists(new_filepath):
		os.mkdir(new_filepath)
	if do_species:
		new_filepath = os.path.join(new_filepath, species_name + '/')
		if not os.path.exists(new_filepath):
			os.mkdir(new_filepath)

	basedata_filepath = os.path.join(DATABASE_PATH, 'introns/' + base_classname + '/')
	if not os.path.exists(basedata_filepath):
		os.mkdir(basedata_filepath)
	if do_species:
		basedata_filepath = os.path.join(basedata_filepath, species_name + '/')
		if not os.path.exists(basedata_filepath):
			os.mkdir(basedata_filepath)
	basedata_filepath = os.path.join(basedata_filepath, 'base_info.dat')
	if not os.path.exists(basedata_filepath):
		if shifted:
			make_shifted_control(base_data)
		else:
			make_shuffled_control(base_data, do_species=do_species, species_name=species_name)
	
	new_base_data = read_base_data(basedata_filepath)

	for ii, data_line in enumerate(base_data):
		bps = int(data_line[0][0])
		start_offset = 0
		end_offset = 0
		if len(data_line[0][1]) > 6:
			start_offset = int(data_line[0][1][-2])
			end_offset = int(data_line[0][1][-1])
		seq = data_line[1]
		new_seq = new_base_data[ii][1]
		new_seq = list(new_seq)
		new_seq[start_offset:(start_offset + 6)] = seq[start_offset:(start_offset + 6)]
		new_seq[(bps-5):(bps+2)] = seq[(bps-5):(bps+2)]
		new_seq[(len(new_seq)-end_offset-3):(len(new_seq)-end_offset)] = \
			seq[(len(new_seq)-end_offset-3):(len(new_seq)-end_offset)]
		new_seq = ''.join(new_seq)
		new_base_data[ii] = (new_base_data[ii][0], new_seq)

	base_data_path = os.path.join(new_filepath, 'base_info.dat')
	write_base_data(base_data_path, new_base_data)

def make_shifted_seq_matched_control(base_data, intron_class):
	base_classname = intron_class + '_shift_' + str(shift_dist)
	make_seq_matched(base_data, base_classname)

def make_shuffle_seq_matched_control(base_data, intron_class, do_species=False, species_name=""):
	base_classname = intron_class + '_shuffle'
	make_seq_matched(base_data, base_classname, shifted=False, \
		do_species=do_species, species_name=species_name)

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Parameters for processing intron data')
	parser.add_argument('intron_class', type=str, help='Intron class that the controls will be based off of')
	parser.add_argument('--make_shuffle', default=False, action='store_true', \
		 help='Make a shuffled intron case-matched control set')
	parser.add_argument('--make_shifted', default=False, action='store_true', \
		 help='Make a shifted intron case-matched control set')
	parser.add_argument('--shift_dist', type=int, default=500, \
		help='Distance to shift the shifted control in the genome')
	parser.add_argument('--make_shifted_seq_matched', default=False, action='store_true', \
		 help='Make a shifted intron case-matched control set with matching splicing sequences')
	parser.add_argument('--make_shuffle_seq_matched', default=False, action='store_true', \
		 help='Make a shuffle intron case-matched control set with matching splicing sequences')
	args = parser.parse_args()

	intron_class = args.intron_class
	make_shuffle = args.make_shuffle
	make_shifted = args.make_shifted
	shift_dist = args.shift_dist
	make_shifted_seq_matched = args.make_shifted_seq_matched
	make_shuffle_seq_matched = args.make_shuffle_seq_matched

	base_data_path = 'introns/' + intron_class + '/base_info.dat'
	base_data_path = os.path.join(DATABASE_PATH, base_data_path)

	base_data = read_base_data(base_data_path)

	if make_shifted:
		make_shifted_control(base_data, intron_class)
	if make_shuffle:
		make_shuffled_control(base_data, intron_class)
	if make_shifted_seq_matched:
		make_shifted_seq_matched_control(base_data, intron_class)
	if make_shuffle_seq_matched:
		make_shuffle_seq_matched_control(base_data, intron_class)
