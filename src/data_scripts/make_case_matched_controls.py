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

def make_shifted_control(base_data):
	new_classname = intron_class + '_shift_' + str(shift_dist)
	new_filepath = os.path.join(DATABASE_PATH, 'introns/' + new_classname + '/')
	if not os.path.exists(new_filepath):
		os.mkdir(new_filepath)

	bed_data = []
	bps = []
	for data_line in base_data: 
		bps += [data_line[0][0]]
		(chr_num, bed_start, bed_end, tag, strand_dir) = data_line[0][1]
		bed_start = int(bed_start)
		bed_end = int(bed_end)
		if strand_dir == '-':
			bed_start = bed_start - shift_dist
			bed_end = bed_end - shift_dist
		else:
			bed_start = bed_start + shift_dist
			bed_end = bed_end + shift_dist
		bed_data += [(chr_num, bed_start, bed_end, tag, strand_dir)]
		
	bed_filename = write_bed(bed_data)

	base_data_path = os.path.join(new_filepath, 'base_info.dat')
	write_base_data_from_bed(base_data_path, bed_filename, bps)
	clear_tmp_files()

def make_shuffled_control(base_data):
	new_classname = intron_class + '_shuffle'
	new_filepath = os.path.join(DATABASE_PATH, 'introns/' + new_classname + '/')
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

def make_seq_matched(base_data, base_classname, shifted=True):
	new_classname = base_classname + '_seq_matched'
	new_filepath = os.path.join(DATABASE_PATH, 'introns/' + new_classname + '/')
	if not os.path.exists(new_filepath):
		os.mkdir(new_filepath)

	basedata_filepath = os.path.join(DATABASE_PATH, 'introns/' + base_classname + '/')
	basedata_filepath = os.path.join(basedata_filepath, 'base_info.dat')
	if not os.path.exists(basedata_filepath):
		if shifted:
			make_shifted_control(base_data)
		else:
			make_shuffled_control(base_data)
	
	new_base_data = read_base_data(basedata_filepath)

	for ii, data_line in enumerate(base_data):
		bps = int(data_line[0][0])
		seq = data_line[1]
		new_seq = new_base_data[ii][1]
		new_seq = list(new_seq)
		new_seq[0:6] = seq[0:6]
		new_seq[(bps-5):(bps+2)] = seq[(bps-5):(bps+2)]
		new_seq[-3:] = seq[-3:]
		new_seq = ''.join(new_seq)
		new_base_data[ii] = (new_base_data[ii][0], new_seq)

	base_data_path = os.path.join(new_filepath, 'base_info.dat')
	write_base_data(base_data_path, new_base_data)

def make_shifted_seq_matched(base_data):
	base_classname = intron_class + '_shift_' + str(shift_dist)
	make_seq_matched(base_data, base_classname)

def make_shuffle_seq_matched(base_data):
	base_classname = intron_class + '_shuffle'
	make_seq_matched(base_data, base_classname, shifted=False)

base_data_path = 'introns/' + intron_class + '/base_info.dat'
base_data_path = os.path.join(DATABASE_PATH, base_data_path)

base_data = read_base_data(base_data_path)

if make_shifted:
	make_shifted_control(base_data)
if make_shuffle:
	make_shuffled_control(base_data)
if make_shifted_seq_matched:
	make_shifted_seq_matched(base_data)
if make_shuffle_seq_matched:
	make_shuffle_seq_matched(base_data)
