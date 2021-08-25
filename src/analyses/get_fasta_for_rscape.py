from core.intron import IntronSet
from config import DATABASE_PATH
from util.gene_file_io import write_fasta
import math

all_introns = IntronSet()
intron_seq_file = DATABASE_PATH + 'introns/standard_allsize/base_info.dat' 
all_introns.init_from_files(intron_seq_file, do_fill_gene_seq=True)
all_introns.fill_ensembl_names()

for intron in all_introns.introns:
	if str(intron.ensembl_name) == 'nan':
		continue
	intron_tag = intron.ensembl_name + "_" + intron.chr_pos[0] + "_" + \
		str(intron.chr_pos[1]) + "_" + str(intron.chr_pos[2])
	filename = DATABASE_PATH + 'introns/intron_flanking/' + intron_tag + '.fa'
	(idx_start, idx_end) = intron.intron_pos_in_gene
	tag = intron.ensembl_name + " " + intron.chr_pos[0] + " " + \
		str(intron.chr_pos[1]) + "-" + str(intron.chr_pos[2]) + " " + \
		str(idx_start) + " " + str(idx_end) + " Saccharomyces cerevisiae S288C"
	print(tag)
	seq = intron.gene_seq

	write_fasta([(tag, seq)], fasta_filename=filename)
