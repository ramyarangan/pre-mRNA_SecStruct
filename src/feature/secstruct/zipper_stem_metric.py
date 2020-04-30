from intron_metric import IntronSet
from intron_metric import Intron

import subprocess # For making RNAfold/RNAcofold calls
import os
import scipy.stats as stats

import matplotlib.pyplot as plt

# This code finds stable stems with dG < DG_LIM that fall within a region that is reasonable
# in the context of the spliceosome. The viable region must be pre-specified based 
# on structural analysis.

# Approach: Starts with reasonably long stems present in the MFE, and then checks the 
# stability of these stems and a slightly expanded surrounding context with RNAfold or RNAcofold. 

# Requires an intron sequence file and secondary structure file. 

# Region is (lowest start, highest end, total distance max, total distance min)
# Stem must start after lowest start + start_idx and must end before end_idx - highest end
DOFIRST=True
first_step_region = (10, 20, 12+10+20, 45+10+20) # For 5'SS to BP
second_step_region = (8, 3, 8+3, 3000) # For BP to 3'SS
species_names = ['agos', 'cgla', 'ecym', 'kafr', 'klac', 'knag', 'kpol', 'kthe', 'kwal', 'ncas', \
	'ndai', 'scer', 'sklu', 'skud', 'smik', 'suva', 'tbla', 'tdel', 'tpha', 'zrou']
region = second_step_region 
if DOFIRST:
	region = first_step_region

# DG_LIM chosen to ensure that at 1 nM two duplexes will hybridize ~50% of the time
DG_LIM=0 # -12 # For a full histogram choose dG = 0
# Minimum intron size
SIZE_MIN=0
# Minimum seed stem size
MIN_SEED=2 # 8
# Prints debug output
DEBUG=False
DEBUG_STEMS=True

# Gets dictionary of start of BP: end of BP
def get_base_pairs(secstruct):
	ii = 0
	stack = []
	bps = {}
	while (ii < len(secstruct)):
		if (secstruct[ii] == "("):
			stack += [ii]
		if (secstruct[ii] == ")"):
			bps[stack[len(stack)-1]] = ii
			stack = stack[:(len(stack)-1)]
		ii += 1
	return bps

# Get stems of the form [(a, b, c, d, n)] for a---b paired to c---d 
# with n bps present in the stem, with max bulge size specified.
# Can handle strings with & in the middle to indicate two strands
def get_stems(secstruct, bulge_size=3):
	bps = get_base_pairs(secstruct)

	stems = []
	ii = 0
	# Search for the start of a new stem
	while (ii < len(secstruct)) and (secstruct[ii] != '&'):
		if ii not in bps:
			ii += 1
			continue
		start = ii
		end = bps[ii]
		jj = bps[ii]
		ii += 1
		bulge_cnt = 0
		n_bps = 1
		# Extend stem as far as possible, as long as no bulge is > limit
		while (ii < jj) and (bulge_cnt < bulge_size + 1):
			# End if you encounter a )
			if secstruct[ii] == ')' or secstruct[ii] == '&':
				break
			# Increment bulge size if necessary
			if secstruct[ii] == '.':
				bulge_cnt += 1
			# If a match is found, restart bulge count and increase bp count
			elif ((bps[ii] <= jj) and (bps[ii] >= jj - bulge_size)):
				jj = bps[ii]
				n_bps += 1
				bulge_cnt = 0
			else: # The endpoint is a different stem
				break

			# Next position
			ii += 1
		stems += [(start, ii - bulge_cnt, jj, end + 1, n_bps)]
	if (DEBUG):
		print(secstruct)
		print("All stems: ")
		print(stems)
	return stems

# Gets location of longest stem with bulge size limit of 3 within this secstruct
def get_max_stem(secstruct, bulge_size=3):
	stems = get_stems(secstruct, bulge_size=bulge_size)
	max_stem = (-1, -1, -1, -1, -1)
	for stem in stems:
		if (stem[4] > max_stem[4]):
			max_stem = stem
	return max_stem

# From the stem locations in 'stems', find all that are compatible with the 
# region limits specified in 'region'
# 'Region' is (lowest start, highest end, total distance max, total distance min)
# Stem must start after lowest start + start_idx and must end before end_idx - highest end
def get_stems_region(stems, region, start_idx, end_idx, seed_stem=MIN_SEED):
	matching_stems = []
	(low, high, max_dist, min_dist) = region
	for stem in stems:
		(a, b, c, d, n_bps) = stem
		if n_bps > seed_stem and (a - start_idx) > low and \
			(end_idx - d) > high and \
			((a - start_idx) + (end_idx - d)) > max_dist and \
			((a - start_idx) + (end_idx - d)) < min_dist:
			matching_stems += [stem]
	if (DEBUG):
		print("Stems in region:")
		print(matching_stems)
	return matching_stems

def collect_dG_secstruct(seq, sys_command):
	f = open('tmp.dat', 'w')
	f.write(seq)
	f.close()
	p = subprocess.Popen(sys_command + ' tmp.dat', shell=True, stdout=subprocess.PIPE)
	lines = p.stdout.readlines()
	os.remove('tmp.dat')

	# String parsing to process RNAfold/cofold output
	secstruct = lines[1].decode("utf-8").split()[0]
	# Handles cases when the output is like ((((..(((((...&)))))..)))) ( -8.20)\n or like
	# ((((..(((((...&)))))..)))) (-8.20)\n
	dG_str = ''.join(lines[1].decode("utf-8").split()[1:]) 
	dG = float(dG_str[1:-1])

	return [dG, secstruct]

# System call to run RNAfold
def run_rnafold(seq):
	return collect_dG_secstruct(seq, 'RNAfold')

# System call to run RNAcofold
def run_cofold(seq):
	return collect_dG_secstruct(seq, 'RNAcofold')

# For an expanded region around the stem candidate from the MFE, check the stability
# using RNAfold or RNAcofold
def get_max_stem_dG(stem, seq, min_start, max_end, expand_low=5, expand_total=25):
	
	# Get region surrounding the stem that does not go too close to the 5'SS and BP
	start1 = max(min_start, stem[0] - expand_low)
	end1 = start1 + expand_total
	end2 = min(max_end, stem[3] + expand_low)
	start2 = end2 - expand_total

	# Assemble input for RNAfold or RNAcofold
	run_seq = ""
	secstruct = ""
	dG = 0
	if start2 > end1 + 4: 	# Two segments are separated in the intron so use RNAcofold
		seq1 = seq[start1:end1]
		seq2 = seq[start2:end2]
		run_seq = seq1 + "&" + seq2
		[dG, secstruct] = run_cofold(run_seq)
	else: # Two segments overlap or near-overlap in the intron; treat as one strand
		run_seq = seq[start1:end2]
		[dG, secstruct] = run_rnafold(run_seq)

	# Get max stem from these segments and make two separate strands
	max_stem = get_max_stem(secstruct)

	# Get dG with RNAcofold
	if (dG < 0):
		seq1 = run_seq[max_stem[0]:max_stem[1]]
		seq2 = run_seq[max_stem[2]:max_stem[3]]
		run_seq = seq1 + "&" + seq2
		[dG, secstruct] = run_cofold(run_seq)

	stem_str = run_seq + "\n" + secstruct + "\n" + str(dG)
	return [dG, stem_str]

def has_stem_dG(bp, seq, mfe, dG_lim=DG_LIM):
	stems = get_stems(mfe)
	stems_region = []
	if not DOFIRST:
		stems_region = get_stems_region(stems, region, bp, len(seq)) # Second step
	else:
		stems_region = get_stems_region(stems, region, 0, bp) # First step
	best_dG = 200 # Some large number
	best_stem = ""
	for stem in stems_region:
		dG = 200
		cur_stem = ""
		if not DOFIRST: 
			[dG, cur_stem] = get_max_stem_dG(stem, seq, bp + region[0], len(seq)-region[1]) # Second step
		else:
			[dG, cur_stem] = get_max_stem_dG(stem, seq, region[0], bp-region[1]) # First step
		if (dG < best_dG):
			best_dG = dG
			best_stem = cur_stem
	return [(best_dG < dG_lim), best_stem, best_dG]

def is_noncanonical(seq):
	is_noncanonical = False
	if seq[0:6] != "GTATGT":
		is_noncanonical = True
	return is_noncanonical

def print_stem_info(seq_filename, secstruct_filename, print_full=False, cutoff=0, dG_lim=-8):
	all_introns = IntronSet()
	all_introns.init_from_files(seq_filename, mfe_filename=secstruct_filename)

	introns_with_stem = []
	num_noncanonical = 0
	num_noncanonical_stem = 0
	num_total = 0
	num_big = 0
	num_big_stem = 0
	for intron in all_introns.introns:
		if len(intron.seq) < cutoff:
			continue
		if len(intron.seq) > 200:
			num_big += 1
		if len(intron.seq) < SIZE_MIN:
			continue
		num_total += 1

		[has_stem, best_stem, best_dG] = has_stem_dG(intron.bp, intron.seq, intron.mfe, dG_lim=dG_lim)

		noncanonical = is_noncanonical(intron.seq)

		if noncanonical:
			num_noncanonical += 1

		if (has_stem):
			if print_full:
				if len(intron.seq) < 400:
					print(intron.name)
					print(intron.seq)
					print(len(intron.seq))
					print(best_stem)
					#print('\n')
			if len(intron.seq) > 200:
				num_big_stem += 1
			introns_with_stem += [intron.name]
			if noncanonical:
				num_noncanonical_stem += 1
	num_stem = len(introns_with_stem)
	num_total = num_total
	summary_str = str(num_stem) + " of " + str(num_total)
	big_str = str(num_big_stem) + " of " + str(num_big)
	noncanonical_str = str(num_noncanonical_stem) + " of " + str(num_noncanonical)
	print(summary_str)
	if print_full:
		_, pval = stats.fisher_exact([[num_noncanonical, num_total - num_noncanonical], \
			[num_noncanonical_stem, num_stem - num_noncanonical_stem]])
		print(noncanonical_str)
		print(big_str)
		print(pval)

def get_dGs_names(seq_filename, secstruct_filename, cutoff=-10, max_num=220):
	all_introns = IntronSet()
	all_introns.init_from_files(seq_filename, mfe_filename=secstruct_filename)
	print(len(all_introns.introns))
	intron_dGs = []
	intron_names = []
	all_names = []
	for ii, intron in enumerate(all_introns.introns):
		if ii > max_num:
			break
		if len(intron.seq) < SIZE_MIN:
			continue
		[has_stem, best_stem, best_dG] = has_stem_dG(intron.bp, intron.seq, intron.mfe, dG_lim=cutoff)

		if (has_stem):
			intron_dGs += [best_dG]
			intron_names += [intron.name]
		all_names += [intron.name]
	return (intron_dGs, intron_names, all_names)

def get_dGs(seq_filename, secstruct_filename, max_num=10000):
	all_introns = IntronSet()
	all_introns.init_from_files(seq_filename, mfe_filename=secstruct_filename)
	print(len(all_introns.introns))
	intron_dGs = []
	for ii, intron in enumerate(all_introns.introns):
		if ii > max_num:
			break
		if len(intron.seq) < SIZE_MIN:
			continue
		[has_stem, best_stem, best_dG] = has_stem_dG(intron.bp, intron.seq, intron.mfe)

		if (has_stem):
			intron_dGs += [best_dG]
	return intron_dGs

def plot_hist_stems(seq_filename, secstruct_filename, control_filename, control_secstruct_filename, \
	plot_lines=[], plot_ext='', max_num=10000):
	intron_dGs = get_dGs(seq_filename, secstruct_filename, max_num=max_num)
	control_dGs = get_dGs(control_filename, control_secstruct_filename, max_num=max_num)
	intron_dGs.sort()
	control_dGs.sort()
	plt.hist([intron_dGs, control_dGs], color=['blue', 'forestgreen'], \
		alpha=0.7, label=['Intron', 'Control'])
	plt.legend(loc='upper left')
	if (len(intron_dGs) > 10):
		plt.axvline(intron_dGs[10], color='blue', linestyle='dashed')
	for plot_line in plot_lines:
		plt.axvline(plot_line, color='black', linewidth=1.0)
	if (len(control_dGs) > 10):
		plt.axvline(control_dGs[10], color='forestgreen', linestyle='dashed')
	name = "BP to 3'SS " + plot_ext
	if DOFIRST:
		name = "5'SS to BP " + plot_ext
	plt.title(name)
	plt.show()

def print_stem_intersects(seq_folder, secstruct_folder):
	scer_seq_filename = seq_folder + 'scer_filtered.dat'
	scer_secstruct_filename = secstruct_folder + 'scer_Vienna_mfe.dat'
	[_, scer_dG_introns, scer_all_introns] = \
		get_dGs_names(scer_seq_filename, scer_secstruct_filename, cutoff=-5)

	for species in species_names:
		seq_filename = seq_folder + species + '_filtered.dat'
		secstruct_filename = secstruct_folder + species + '_Vienna_mfe.dat'
		[dGs, dG_introns, all_introns] = get_dGs_names(seq_filename, secstruct_filename, cutoff=-5)
		not_scer = 0 # Introns in these species with stems where the intron is lost in scer
		not_stem_scer = 0 # Introns in these species without a stem in scer but with intron present
		total_passing_cutoff = 0
		for ii, dG_intron in enumerate(dG_introns):
			if dGs[ii] > -10:
				continue
			total_passing_cutoff += 1
			in_scer = False
			for intron in scer_all_introns:
				if intron == dG_intron:
					in_scer = True
			if not in_scer:
				not_scer += 1
			if in_scer:
				stem_in_scer = False
				for intron in scer_dG_introns:
					if intron == dG_intron:
						stem_in_scer = True
				if not stem_in_scer:
					not_stem_scer += 1
		print(species)
		print('Total number of stem-containing introns: ' + str(total_passing_cutoff))
		print('Not a intron-containing gene in scer: ' + str(not_scer))
		print('Not a stem-containing intron in scer: ' + str(not_stem_scer))


def plot_hist_stems_species(seq_folder, control_folder, secstruct_folder):
	for species in species_names:
		plot_hist_stems(seq_folder + species + '_filtered.dat', \
			secstruct_folder + species + '_Vienna_mfe.dat', \
			control_folder + species + '_control_filtered.dat', \
			secstruct_folder + species + '_control_Vienna_mfe.dat', \
			plot_ext=species)

def plot_hist_stems_all(seq_filename, secstruct_filename, control_filename1, control_secstruct_filename1, 
		control_filename2, control_secstruct_filename2, control_filename3, control_secstruct_filename3, max_num=220):
	intron_dGs = get_dGs(seq_filename, secstruct_filename, max_num=max_num)
	control_dGs1 = get_dGs(control_filename1, control_secstruct_filename1, max_num=max_num)
	control_dGs2 = get_dGs(control_filename2, control_secstruct_filename2, max_num=max_num)
	control_dGs3 = get_dGs(control_filename3, control_secstruct_filename3, max_num=max_num)
	intron_dGs.sort()
	control_dGs1.sort()
	control_dGs2.sort()
	control_dGs3.sort()
	print(len(control_dGs1))
	print(len(control_dGs2))
	plt.hist([intron_dGs, control_dGs1, control_dGs2, control_dGs3], color=['blue', 'orange', 'red', 'forestgreen'], \
		alpha=0.7, label=['Intron', 'Fake Standard Introns', 'Fake Proto Introns', 'Control'])
	plt.legend(loc='upper left')
	if (len(intron_dGs) > 10):
		plt.axvline(intron_dGs[10], color='blue', linestyle='dashed')
	if (len(control_dGs1) > 10):
		plt.axvline(control_dGs1[10], color='orange', linestyle='dashed')
	if (len(control_dGs1) > 10):
		plt.axvline(control_dGs2[10], color='red', linestyle='dashed')
	if (len(control_dGs2) > 10):
		plt.axvline(control_dGs3[10], color='forestgreen', linestyle='dashed')
	name = "BP to 3'SS"
	if DOFIRST:
		name = "5'SS to BP"
	plt.title(name)
	plt.show()

def print_stems_scer():
	# print_stem_info('../data/seqs/all_introns_filtered.dat', 
	# 	'../data/secstructs/introns/Vienna_mfe.dat', 
	# 	print_full=True)
	# print_stem_info('../data/seqs/control_500_filtered.dat', 
	# 	'../data/secstructs/control_500/Vienna_mfe.dat')
	print_stem_info('../data/seqs/ares/standard_introns_filtered.dat', 
		'../data/secstructs/ares/introns/Vienna_mfe.dat', 
		print_full=True)
	print_stem_info('../data/seqs/ares/standard_control_500_filtered.dat', 
		'../data/secstructs/ares/control_500/Vienna_mfe.dat')
	print_stem_info('../data/seqs/ares/standard_control_500_matchseq_filtered.dat', 
		'../data/secstructs/ares/control_500_matchseq/Vienna_mfe.dat')
	print_stem_info('../data/seqs/ares/standard_control_shuffled_filtered.dat', 
		'../data/secstructs/ares/control_shuffled/Vienna_mfe.dat')

def print_stems_scer_kthe():
	print_stem_info('../data/seqs/species/seqs_matched/scer_filtered.dat', 
		'../data/secstructs/species/scer_Vienna_mfe.dat', 
		print_full=True)
	print_stem_info('../data/seqs/species/seqs_matched/kthe_filtered.dat', 
		'../data/secstructs/species/kthe_Vienna_mfe.dat', 
		print_full=True)

def plot_hist_proto_standard_decoy():
	plot_hist_stems_all('../secstructs/data/seqs/ares/standard_introns_filtered.dat', 
		'../secstructs/data/secstructs/ares/introns/Vienna_mfe.dat', 
		'fake_standard_filtered.dat', 'standard_Vienna_mfe.dat', 
		'fake_proto_filtered.dat', 'proto_Vienna_mfe.dat', 
		'../secstructs/data/seqs/ares/standard_control_500_filtered.dat', 
		'../secstructs/data/secstructs/ares/control_500/Vienna_mfe.dat')

def do_pombe_analysis():
	print_stem_info('../data/seqs/species/pombe/pombe_filtered.dat',
		'../data/secstructs/species/pombe_Vienna_mfe.dat', print_full=True, cutoff=100, dG_lim=-10)
	print_stem_info('../data/seqs/species/pombe/pombe_control_filtered.dat',
		'../data/secstructs/species/pombe_control_Vienna_mfe.dat', print_full=False, cutoff=100, dG_lim=-10)
	plot_hist_stems('../data/seqs/species/pombe/pombe_filtered.dat',
		'../data/secstructs/species/pombe_Vienna_mfe.dat',
		'../data/seqs/species/pombe/pombe_control_filtered.dat',
		'../data/secstructs/species/pombe_control_Vienna_mfe.dat')

def do_species_analysis():
	plot_hist_stems_species('../data/seqs/species/seqs_matched/', \
		'../data/seqs/species/control_matched/',\
		'../data/secstructs/species/')
	print_stem_intersects('../data/seqs/species/seqs_matched/',	'../data/secstructs/species/')

def plot_scer_stem_hist():
	plot_hist_stems('../data/seqs/ares/standard_introns_filtered.dat', 
	 	'../data/secstructs/ares/introns/Vienna_mfe.dat', 
	 	'../data/seqs/ares/standard_control_500_filtered.dat', 
	 	'../data/secstructs/ares/control_500/Vienna_mfe.dat', 
	 	plot_lines=[-15.3, -12.2, -14.7, -15.6, -14.6, -17.6, -15.9, \
	 	-12.6, -12.8, -12.7, -17.4, -13.1])

print_stems_scer()
# do_pombe_analysis()
