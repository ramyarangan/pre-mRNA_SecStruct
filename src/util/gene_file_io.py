import random 
import os
from config import TMP_PATH
from config import DATABASE_PATH

def get_rnd_tmp_file():
	return TMP_PATH + str(int(random.random() * 10000))

def clear_tmp_files(): 
	for f in os.listdir(TMP_PATH):
		os.remove(os.path.join(TMP_PATH, f))

# Expects list of tuples with bed file info: 
# (chr_num, bed_start, bed_end, tag, strand_dir)
def write_bed(bed_data):
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
			bed_line.split('\t')
		bed_items += [(chr_num, bed_start, bed_end, tag, strand_dir)]
	return bed_items

def read_base_data(basedata_filename):
	f = open(basedata_filename)
	all_lines = f.readlines()
	f.close()

	# A list of tuples of the form: 
	# ((bp, (chr_num, bed_start, bed_end, tag, strand_dir)), sequence)
	base_data = []
	for ii in range(int(len(all_lines)/2)):
		bp = all_lines[ii * 2].split(' ')[0]
		bed_items = all_lines[ii * 2].split(' ')[1].split('\t')
		bed_items = tuple(bed_items[0:4] + bed_items[5:6])
		base_data += [((bp, bed_items), seq)]
	
	return base_data

# A base data file can be written using the bed file and branchpoint positions
def write_base_data_from_bed(base_data_path, bed_filename, bps):
	genome_file = DATABASE_PATH + 'genome/sacCer3.fa'
	fasta_filename = fasta_from_bed(bedfile, genome_file)
	fasta_seqs = read_fasta(fasta_filename)
	bed_items = read_bedfile(bed_filename)
	os.remove(fasta_filename)

	if len(bed_items) != len(fasta_seqs):
		raise RuntimeError("Fasta sequence not recovered for some bed lines")

	f = open(base_data_path, 'w')
	for ii, (_, seq) in enumerate(fasta_seqs): 
		f.write("%d %s\n" % (int(bps[ii]), '\t'.join(bed_items[1:])))
		f.write("%s\n" % seq)
	f.close()

# base_info_items: A list of tuples of the form: 
# (branchpoint position, chr name, chr start idx, \
# chr end idx, gene name, 0, strand_dir)
def write_base_data(base_info_path, base_info_items, base_info_seqs):
	f = open(base_info_path, 'w')
	for ii, base_info_item in enumerate(base_info_items):
		f.write("%d %s\n" % (base_info_item[0], '\t'.join(base_info_item[1:])))
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

