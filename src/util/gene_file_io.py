import random 
import os
from config import TMP_PATH

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




