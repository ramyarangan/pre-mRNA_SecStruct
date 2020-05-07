class Intron:
	def __init__(self, seq='', bp=-1, mfe='', ens=[], name='', chr_pos=('',-1,-1), strand='-'):
		self.bp = bp
		self.seq = seq
		self.mfe = mfe
		self.ens = ens
		self.name = name
		self.chr_pos = chr_pos
		self.strand = strand

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
			seq = seq_lines[2 * ii + 1].split()[0]
			mfe = ''
			ens = []
			if (len(mfe_lines) != 0):
				mfe = mfe_lines[2 * ii + 1].split()[0]
			if (len(ens_lines) != 0):
				ens = ens_lines[((ens_size+1) * ii + 1):((ens_size+1) * (ii + 1))]
				ens = [x.split()[0] for x in ens]

			self.introns.append(Intron(seq, bp, mfe, ens, name, chr_pos, strand))
