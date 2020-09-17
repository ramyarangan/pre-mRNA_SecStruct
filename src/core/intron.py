from util.gene_names import get_ensembl_names

class Intron:
	def __init__(self, seq='', bp=-1, mfe='', ens=[], name='', 
		chr_pos=('',-1,-1), strand='-', ensembl_name='', 
		fivess_offset=0, threess_offset=0):
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

class IntronSet:
	def __init__(self):
		self.introns = []

	def init_from_files(self, seq_filename, mfe_filename="", \
					ens_filename="", ens_size=1000):
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
			intron_info = seq_lines[2 * ii].split()
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
				chr_pos, strand, fivess_offset, threess_offset))

		self.fill_ensembl_names()

	def fill_ensembl_names(self):
		refseq_names = [intron.name for intron in self.introns]
		ensembl_genes = get_ensembl_names(refseq_names)
		for ii, intron in enumerate(self.introns):
			intron.ensembl_name = ensembl_genes[ii]

	def get_intron_dict(self):
		intron_dict = {}
		for intron in self.introns:
			intron_dict[intron.name] = intron
		return intron_dict

