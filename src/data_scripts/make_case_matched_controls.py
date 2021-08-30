# 8/27/21
# Creates shuffled and shifted controls given an intron class
# Each intron in the base_info.dat file is line-matched to the original
# intron class (order matters!)

import os
import argparse
import random 
from config import DATABASE_PATH
from util.gene_file_io import * 

parser.add_argument('intron_class', type=str, help='Intron class that the controls will be based off of')
parser.add_argument('--make_shuffle', default=False, action='store_true', \
	 help='Make a shuffled intron case-matched control set')
parser.add_argument('--make_shifted', default=False, action='store_true', \
	 help='Make a shifted intron case-matched control set')
parser.add_argument('--shift_dist', help='Make a shifted intron case-matched control set', 
	type=int, default=500, help='Distance to shift the shifted control in the genome')
args = parser.parse_args()

intron_class = args.intron_class
make_shuffle = args.make_shuffle
make_shifted = args.make_shifted
shift_dist = args.shift_dist

def make_shifted_control(base_data):
	new_classname = intron_class + '_shift_' + str(shift_dist)
	new_filepath = os.path.join(DATABASE_PATH, 'introns/' + new_classname + '/')
	os.mkdir(new_filepath)

	bed_data = []
	bps = []
	for data_line in base_data: 
		bps += [data_line[0][0]]
		(chr_num, bed_start, bed_end, tag, strand_dir) = data_line[0][1]
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
	os.mkdir(new_filepath)
	base_data_path = os.path.join(new_filepath, 'base_info.dat')

	base_data_items = []
	base_data_seqs = []
	for data_line in base_data: 
		((bp, bed_items), seq) = data_line
		base_data_items += [tuple(list(bp) + list(bed_items))]
		base_data_seqs += ''.join(random.sample(seq, len(seq)))
	
	write_base_data(base_data_path, base_data_items, base_data_seqs)

base_data_path = 'introns/' + intron_class + '/base_info.dat'
base_data_path = os.path.join(DATABASE_PATH, base_data_path)

base_data = read_base_data(base_data_path)

if make_shifted:
	make_shifted_control(base_data)
if make_shuffled:
	make_shuffled_control(base_data)
