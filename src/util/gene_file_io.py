import random 
import os
from config import TMP_PATH
from config import DATABASE_PATH
import numpy as np 

def get_rnd_tmp_file():
	return TMP_PATH + str(int(random.random() * 10000))

def clear_tmp_files(): 
	for f in os.listdir(TMP_PATH):
		os.remove(os.path.join(TMP_PATH, f))

# Expects list of tuples with bed file info: 
# (chr_num, bed_start, bed_end, tag, strand_dir)
def write_bed(bed_data, filename=""):
	bed_filename = filename
	if filename == "":
		bed_filename = get_rnd_tmp_file() + ".bed"
	
	f = open(bed_filename, 'w')
	for (chr_num, bed_start, bed_end, tag, strand_dir) in bed_data: 
		bedline = '\t'.join([chr_num, str(bed_start), str(bed_end), \
			tag, "0", strand_dir])
		f.write("%s\n" % bedline)
	f.close()

	return bed_filename

def read_bed(bed_filename):
	f = open(bed_filename)
	bed_lines = f.readlines()
	f.close()

	bed_items = []
	for bed_line in bed_lines:
		(chr_num, bed_start, bed_end, tag, _, strand_dir) = \
			bed_line.replace('\n', '').split('\t')
		bed_items += [(chr_num, bed_start, bed_end, tag, strand_dir)]
	return bed_items

def read_base_data(basedata_filename):
	f = open(basedata_filename)
	all_lines = f.readlines()
	f.close()

	# A list of tuples of the form: 
	# ((bp, (chr_num, bed_start, bed_end, tag, strand_dir)), sequence)
	# or
	# ((bp, (chr_num, bed_start, bed_end, tag, strand_dir, fivess_offset, threess_offset)), sequence)
	base_data = []
	for ii in range(int(len(all_lines)/2)):
		bp = all_lines[ii * 2].split(' ')[0]
		bed_items = all_lines[ii * 2].split(' ')[1].strip('\n').split('\t')
		bed_items = tuple(bed_items[0:4] + bed_items[5:])
		seq = all_lines[ii * 2 + 1].replace('\n', '')
		base_data += [((bp, bed_items), seq)]
	
	return base_data

def write_base_data(base_data_path, base_data):
	f = open(base_data_path, 'w')
	for base_item in base_data:
		bed_first = '\t'.join(base_item[0][1][0:4])
		bed_second = '\t'.join(base_item[0][1][4:])
		f.write("%d %s\t0\t%s\n" % (int(base_item[0][0]), bed_first, bed_second))
		f.write("%s\n" % base_item[1])
	f.close()

# A base data file can be written using the bed file and branchpoint positions
def write_base_data_from_bed(base_data_path, bed_filename, bps, offsets=[]):
	genome_file = DATABASE_PATH + 'genome/sacCer3.fa'
	fasta_filename = fasta_from_bed(bed_filename, genome_file)
	fasta_seqs = read_fasta(fasta_filename)
	bed_items = read_bed(bed_filename)
	os.remove(fasta_filename)

	if len(bed_items) != len(fasta_seqs):
		raise RuntimeError("Fasta sequence not recovered for some bed lines")

	f = open(base_data_path, 'w')
	for ii, (_, seq) in enumerate(fasta_seqs):
		bed_first = '\t'.join(bed_items[ii][:-1])
		bed_second = bed_items[ii][-1]
		if len(offsets) > 0:
			bed_second = '\t'.join([bed_second] + list(offsets[ii]))
		f.write("%d %s\t0\t%s\n" % (int(bps[ii]), bed_first, bed_second))
		f.write("%s\n" % seq)
	f.close()

# base_info_items: A list of tuples of the form: 
# (branchpoint position, chr name, chr start idx, \
# chr end idx, gene name, 0, strand_dir)
def write_base_data_items(base_info_path, base_info_items, base_info_seqs):
	f = open(base_info_path, 'w')
	for ii, base_info_item in enumerate(base_info_items):
		base_info_item_to_str = [str(x) for x in base_info_item[1:]]
		bed_first = '\t'.join(base_info_item_to_str[:-1])
		bed_second = base_info_item_to_str[-1]
		f.write("%d %s\t0\t%s\n" % (int(base_info_item[0]), bed_first, bed_second))
		f.write("%s\n" % base_info_seqs[ii])
	f.close()

def fasta_from_bed(bedfile, genome_fasta):
	fasta_filename = get_rnd_tmp_file() + ".fa"

	flags = '-fi ' + genome_fasta + ' -bed ' + bedfile + \
		' -s -name -fo ' + fasta_filename
	os.system('fastaFromBed ' + flags)
	return fasta_filename

def read_fasta(fasta_filename):
	f = open(fasta_filename)
	all_lines = f.readlines()
	f.close()

	cur_tag = None
	cur_seq = ""
	fasta_seqs = []
	for line in all_lines:
		if line[0] == ">":
			if cur_tag is not None:
				fasta_seqs += [(cur_tag, cur_seq)]
			cur_tag = line[1:]
			cur_seq = ""
		else:
			cur_seq += line.replace('\n', '')
	if cur_tag is not None:
		fasta_seqs += [(cur_tag, cur_seq)]

	return fasta_seqs

# Expects list of tuples with: 
# (tag, sequence)
def write_fasta(fasta_seqs, fasta_filename="", line_break=50):
	if fasta_filename == "":
		fasta_filename = get_rnd_tmp_file() + ".fa"

	f = open(fasta_filename, 'w')
	for (tag, seq) in fasta_seqs:
		tag_line = ">" + tag
		f.write("%s\n" % tag_line)

		for ii in range(1 + int(len(seq)/line_break)):
			end_idx = min((ii + 1) * line_break, len(seq))
			fasta_window = seq[(ii * line_break):end_idx]
			f.write("%s\n" % fasta_window)

	f.close()

def read_zipper_stem_file(zipper_stem_filename):
	f = open(zipper_stem_filename)
	zipper_stem_lines = f.readlines()
	f.close()

	zipper_stem_data = []
	for ii in range(int(len(zipper_stem_lines)/3)):
		seq = zipper_stem_lines[ii * 3].replace('\n', '')
		seq1 = seq.split("&")[0]
		seq2 = seq.split("&")[1]
		secstruct = zipper_stem_lines[ii * 3 + 1].replace('\n', '')
		dG = zipper_stem_lines[ii * 3 + 2].replace('\n', '')

		zipper_stem_data += [(seq1, seq2, secstruct, dG)]
		
	return zipper_stem_data

def write_zipper_stem_file(zipper_stem_filename, zipper_stem_data):
	f = open(zipper_stem_filename, 'w')

	for cur_data in zipper_stem_data: 
		[seq1, seq2, secstruct, dG] = cur_data
		f.write("%s&%s\n%s\n%f\n" % (seq1, seq2, secstruct, dG))

	f.close()


def process_zipper_stem_entry(has_dG, best_stem):
	if not has_dG: 
		return ("", "", "", 0)

	seq1 = best_stem.split()[0].split("&")[0]
	seq2 = best_stem.split()[0].split("&")[1]
	secstruct = best_stem.split()[1]
	dG = float(best_stem.split()[2])

	return (seq1, seq2, secstruct, dG)

def get_bpp_from_file(bpp_dir, chr_pos, seq_len):
	(chr_num, start_idx, end_idx) = chr_pos

	filename = chr_num + "_" + str(start_idx) + "_" + str(end_idx) + "_bpp.csv"
	bpp_file = os.path.join(bpp_dir, filename)

	f = open(bpp_file)
	bpp_lines = f.readlines()
	f.close()

	bpp_arr = np.zeros((seq_len, seq_len))
	for ii in range(seq_len):
		bpp_items = bpp_lines[ii].replace('\n', '').split(',')
		bpp_arr[ii,:] = bpp_items

	return bpp_arr

