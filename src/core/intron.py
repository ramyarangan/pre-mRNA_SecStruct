from util.gene_names import get_ensembl_names
from core.gene import GeneSet
from util.general import reverse_invert
from config import DATABASE_PATH
from util.gene_file_io import *

class Intron:
	def __init__(self, seq='', bp=-1, mfe='', ens=[], name='', 
		chr_pos=('',-1,-1), strand='-', ensembl_name='', 
		fivess_offset=0, threess_offset=0, 
		gene_seq="", gene_start=-1, gene_end=-1, intron_pos_in_gene=-1, 
		do_fill_gene_seq=False):

		self.bp = bp
		self.seq = seq
		self.mfe = mfe
		self.ens = ens
		self.name = name # refseq name
		self.chr_pos = chr_pos
		self.strand = strand
		self.ensembl_name = ensembl_name
		self.fivess_offset = fivess_offset
		self.threess_offset = threess_offset

		# Initialize gene values
		# For inter-gene introns, include the full gene sequence
		# For 5'UTR introns, include 10 nts before the intron and the full gene 
		self.gene_seq = gene_seq 
		self.gene_start = gene_start # Start position of gene
		self.gene_end = gene_end # End position of gene
		self.intron_pos_in_gene = intron_pos_in_gene # (Start, end) position of intron within the gene sequence
		if do_fill_gene_seq: 
			self.fill_gene_seq()

	def fill_gene_seq(self, UTR_offset=10):
		print(self.name)
		gene_set = GeneSet()
		genes_dict = gene_set.get_genes_dict()
		if self.name not in genes_dict:
			# Figure out why some gene annotations aren't right?
			print("Intron not in ORF: %s" % (self.name))
			return

		gene = genes_dict[self.name]
		
		chr_num = gene.chr_num
		(chr_start, chr_end) = gene.chr_pos
		strand_dir = gene.strand_dir
		self.gene_start = chr_start
		self.gene_end = chr_end

		bed_start = chr_start
		bed_end = chr_end

		(_, intron_start, intron_end) = self.chr_pos
		if intron_start < chr_start: 
			bed_start = intron_start - UTR_offset
		if intron_end > chr_end: 
			bed_end = intron_end + UTR_offset

		bedfile = write_bed([(chr_num, bed_start, bed_end, self.name, strand_dir)])
		genome_file = DATABASE_PATH + 'genome/sacCer3.fa'
		fasta_file = fasta_from_bed(bedfile, genome_file)
		fasta_seqs = read_fasta(fasta_file)
		(tag, gene_seq) = fasta_seqs[0]

		if self.seq not in gene_seq:
			print("Gene sequence: " + gene_seq)
			print("Intron sequence: " + self.seq)
			raise RuntimeError("Intron not in gene sequence")

		self.gene_seq = gene_seq
		idx_start = gene_seq.index(self.seq)
		idx_end = idx_start + len(self.seq)
		self.intron_pos_in_gene = (idx_start, idx_end)
		clear_tmp_files()

class IntronSet:
	def __init__(self):
		self.introns = []

	def init_from_files(self, seq_filename, mfe_filename="", \
					ens_filename="", ens_size=1000, do_fill_gene_seq=False):
		f = open(seq_filename)
		seq_lines = f.readlines()
		f.close()

		mfe_lines = []
		if (mfe_filename != ""):
			f = open(mfe_filename)
			mfe_lines = f.readlines()
			f.close()
		ens_lines = []
		if (ens_filename != ""):
			f = open(ens_filename)
			ens_lines = f.readlines()
			f.close()

		self.introns = []
		for ii in range(int(len(seq_lines)/2)):
			intron_info = seq_lines[2 * ii].split(' ')[1].split('\t')
			intron_info = [seq_lines[2 * ii].split(' ')[0]] + intron_info
			intron_info = [x.replace('\n', '') for x in intron_info]
			bp = int(intron_info[0])
			chr_pos = (intron_info[1], int(intron_info[2]), int(intron_info[3]))
			name = intron_info[4]
			strand = intron_info[6]
			fivess_offset = 0
			threess_offset = 0
			if len(intron_info) > 7:
				fivess_offset = int(intron_info[7])
				threess_offset = int(intron_info[8])

			seq = seq_lines[2 * ii + 1].split()[0]
			mfe = ''
			ens = []
			if (len(mfe_lines) != 0):
				mfe = mfe_lines[2 * ii + 1].split()[0]
			if (len(ens_lines) != 0):
				ens = ens_lines[((ens_size+1) * ii + 1):((ens_size+1) * (ii + 1))]
				ens = [x.split()[0] for x in ens]

			self.introns.append(Intron(seq, bp, mfe, ens, name, \
				chr_pos, strand, fivess_offset, threess_offset, \
				do_fill_gene_seq=do_fill_gene_seq))

		self.fill_ensembl_names()

	def fill_ensembl_names(self):
		refseq_names = [intron.name for intron in self.introns]
		no_blank_refseq_names = []
		for refseq_name in refseq_names:
			# Some introns haven't been properly associated with a gene name
			# because they were in the UTR region outside the gene annotated window
			if refseq_name == "":
				continue
			no_blank_refseq_names += [refseq_name]
		ensembl_genes = get_ensembl_names(no_blank_refseq_names)
		cnt = 0
		for ii, intron in enumerate(self.introns):
			if intron.name == "":
				continue
			intron.ensembl_name = ensembl_genes[cnt]
			cnt += 1

	def get_intron_dict(self):
		intron_dict = {}
		for intron in self.introns:
			intron_dict[intron.name] = intron
		return intron_dict

