from core.intron import *
from config import DATABASE_PATH
from arnie.sample_structures import sample_structures
from arnie.mfe import mfe

MFE_PACKAGE_LIST = ['vienna', 'contrafold', 'rnastructure', 'eternafold']
ENS_PACKAGE_LIST = ['vienna', 'rnastructure']

def add_secstruct_mfe_to_database(intron_class, sec_struct_pkg, print_freq=10):
	intron_seq_file = DATABASE_PATH + 'introns/' + intron_class + '/base_info.dat' 
	all_introns = IntronSet()
	all_introns.init_from_files(intron_seq_file, get_ensembl_names=False)

	secstruct_file = DATABASE_PATH + 'introns/' + intron_class + '/' + sec_struct_pkg + '_mfe.dat'

	f = open(secstruct_file, 'w')
	for ii, intron in enumerate(all_introns.introns):
		if ii % print_freq == 0: 
			print("Done with %d/%d intron MFE structure predictions." % \
				(ii, len(all_introns.introns)))

		if sec_struct_pkg.lower() not in MFE_PACKAGE_LIST:
			raise NotImplementedError()
		mfe_struct = mfe(intron.seq.replace('T', 'U'), package=sec_struct_pkg.lower())
		
		bp = intron.bp
		(chr_name, chr_start, chr_end) = intron.chr_pos
		name = intron.name
		strand = intron.strand

		f.write(">%d %s %d %d %s 0 %s\n" % \
			(bp, chr_name, chr_start, chr_end, name, strand))
		f.write("%s\n" % mfe_struct)
	f.close()

def add_secstruct_ens_to_database(intron_class, sec_struct_pkg, print_freq=10):
	intron_seq_file = DATABASE_PATH + 'introns/' + intron_class + '/base_info.dat' 
	all_introns = IntronSet()
	all_introns.init_from_files(intron_seq_file, get_ensembl_names=False)

	secstruct_file = DATABASE_PATH + 'introns/' + intron_class + '/' + sec_struct_pkg + '_ens.dat'

	f = open(secstruct_file, 'w')
	for ii, intron in enumerate(all_introns.introns):
		if ii % print_freq == 0: 
			print("Done with %d/%d intron ensemble structure predictions." % \
				(ii, len(all_introns.introns)))

		if sec_struct_pkg.lower() not in ENS_PACKAGE_LIST:
			raise NotImplementedError()
		ens = sample_structures(intron.seq.replace('T', 'U'), n_samples=1000, \
			package=sec_struct_pkg.lower())
		
		bp = intron.bp
		(chr_name, chr_start, chr_end) = intron.chr_pos
		name = intron.name
		strand = intron.strand

		f.write(">%d %s %d %d %s 0 %s\n" % \
			(bp, chr_name, chr_start, chr_end, name, strand))
		for struct in ens:
			f.write("%s\n" % struct)
	f.close()
