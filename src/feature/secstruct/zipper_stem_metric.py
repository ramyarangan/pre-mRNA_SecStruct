from feature.secstruct.secstruct_metric import SecstructMetric

import subprocess # For making RNAfold/RNAcofold calls
import os
import scipy.stats as stats

import matplotlib.pyplot as plt


# This code finds stable stems with dG < DG_LIM that fall within a region that is reasonable
# in the context of the spliceosome. The viable region must be pre-specified based 
# on structural analysis.

# Approach: Starts with reasonably long stems present in the MFE, and then checks the 
# stability of these stems and a slightly expanded surrounding context with RNAfold or RNAcofold. 

class ZipperStemMetric(SecstructMetric):
	def __init__(self, do_first=True, dg_lim=0, min_seed=0, name="ZipperStemMetric", max_ens=100):
		# Region is (lowest start, highest end, total distance max, total distance min)
		# Stem must start after lowest start + start_idx and must end before end_idx - highest end
		self.DOFIRST=do_first
		first_step_region = (10, 20, 12+10+20, 55+10+20) # For 5'SS to BP
		second_step_region =  (0, 0, 0, 3000)# For BP to 3'SS (8, 3, 8+3, 3000)
		self.region = second_step_region 
		if self.DOFIRST:
			self.region = first_step_region
		self.name = name
		self.DG_LIM=dg_lim# -12 # For a full histogram choose dG = 0
		self.MIN_SEED=min_seed # 8 # Minimum seed stem size
		self.max_ens = max_ens # Maximum number of secondary structures to average over when doing
								# ensemble calculations

	# Gets dictionary of start of BP: end of BP
	def get_base_pairs(self, secstruct):
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
	def get_stems(self, secstruct, bulge_size=5):
		bps = self.get_base_pairs(secstruct)

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
		return stems

	# Gets location of longest stem with bulge size limit of 3 within this secstruct
	def get_max_stem(self, secstruct, bulge_size=3):
		stems = self.get_stems(secstruct, bulge_size=bulge_size)
		max_stem = (-1, -1, -1, -1, -1)
		for stem in stems:
			if (stem[4] > max_stem[4]):
				max_stem = stem
		return max_stem

	# From the stem locations in 'stems', find all that are compatible with the 
	# region limits specified in 'region'
	# 'Region' is (lowest start, highest end, total distance max, total distance min)
	# Stem must start after lowest start + start_idx and must end before end_idx - highest end
	def get_stems_region(self, stems, region, start_idx, end_idx):
		seed_stem = self.MIN_SEED
		matching_stems = []
		(low, high, max_dist, min_dist) = region
		for stem in stems:
			(a, b, c, d, n_bps) = stem
			if n_bps > seed_stem and (a - start_idx) > low and \
				(end_idx - d) > high and \
				((a - start_idx) + (end_idx - d)) > max_dist and \
				((a - start_idx) + (end_idx - d)) < min_dist:
				matching_stems += [stem]
		return matching_stems

	def collect_dG_secstruct(self, seq, sys_command):
		f = open('tmp.dat', 'w')
		f.write(seq.replace('T', 'U'))
		f.close()
		try: 
			p = subprocess.Popen(sys_command + ' tmp.dat', shell=True, stdout=subprocess.PIPE)
			lines = p.stdout.readlines()
			os.remove('tmp.dat')

			# String parsing to process RNAfold/cofold output
			secstruct = lines[1].decode("utf-8").split()[0]
			# Handles cases when the output is like ((((..(((((...&)))))..)))) ( -8.20)\n or like
			# ((((..(((((...&)))))..)))) (-8.20)\n
			dG_str = ''.join(lines[1].decode("utf-8").split()[1:]) 
			dG = float(dG_str[1:-1])
		except: 
			dG = self.DG_LIM
			secstruct = ""
		return [dG, secstruct]

	# System call to run RNAfold
	def run_rnafold(self, seq):
		return self.collect_dG_secstruct(seq, 'RNAfold')

	# System call to run RNAcofold
	def run_cofold(self, seq):
		return self.collect_dG_secstruct(seq, 'RNAcofold')

	# For an expanded region around the stem candidate from the MFE, check the stability
	# using RNAfold or RNAcofold
	def get_max_stem_dG(self, stem, seq, min_start, max_end, expand_low=5): #, expand_total=25):
	
		# Get region surrounding the stem that does not go too close to the 5'SS and BP
		start1 = max(min_start, stem[0] - expand_low)
		end1 = stem[1] + expand_low
		end2 = min(max_end, stem[3] + expand_low)
		start2 = stem[2] - expand_low

		# Assemble input for RNAfold or RNAcofold
		run_seq = ""
		secstruct = ""
		dG = self.DG_LIM
		seq1 = seq[start1:min(end1, start2)]
		seq2 = seq[max(end1, start2):end2]
		run_seq = seq1 + "&" + seq2
		if len(run_seq) > 2:
			[dG, secstruct] = self.run_cofold(run_seq)

		stem_str = run_seq + "\n" + secstruct + "\n" + str(dG)
		return [dG, stem_str]
		
		"""
		Old function definition
		# Get region surrounding the stem that does not go too close to the 5'SS and BP
		start1 = max(min_start, stem[0] - expand_low)
		end1 = start1 + expand_total
		end2 = min(max_end, stem[3] + expand_low)
		start2 = end2 - expand_total

		# Assemble input for RNAfold or RNAcofold
		run_seq = ""
		secstruct = ""
		dG = self.DG_LIM
		if start2 > end1 + 4: 	# Two segments are separated in the intron so use RNAcofold
			seq1 = seq[start1:end1]
			seq2 = seq[start2:end2]
			run_seq = seq1 + "&" + seq2
			if len(run_seq) > 2:
				[dG, secstruct] = self.run_cofold(run_seq)
		else: # Two segments overlap or near-overlap in the intron; treat as one strand
			run_seq = seq[start1:end2]
			if len(run_seq) > 1:
				[dG, secstruct] = self.run_rnafold(run_seq)

		# Get max stem from these segments and make two separate strands
		max_stem = self.get_max_stem(secstruct)

		# Get dG with RNAcofold
		if (dG < self.DG_LIM):
			seq1 = run_seq[max_stem[0]:max_stem[1]]
			seq2 = run_seq[max_stem[2]:max_stem[3]]
			run_seq = seq1 + "&" + seq2
			if len(run_seq) > 2:
				[dG, secstruct] = self.run_cofold(run_seq)

		stem_str = run_seq + "\n" + secstruct + "\n" + str(dG)
		return [dG, stem_str]
		"""

	# Get the highest average BPP for a helix of length at least 4 within the stem
	def get_best_bpp(self, secstruct, stem, bpp_matrix, len_cutoff=4, bpp_thresh=0.7):
		
		strand1 = secstruct[stem[0]:stem[1]]
		strand2 = secstruct[stem[2]:stem[3]]

		bps = get_base_pairs(secstruct)

		ii = stem[0]
		cur_end = stem[3] - 1
		cur_cnt = 0
		total_bpp = 0
		bpps = []
		cnts = []
		while ii < stem[1] and cur_end >= stem[2]:
			if (ii in bps) and (bps[ii] == cur_end):
				cur_cnt += 1
				# For now get the maximum base-pairing probability in the stem 
				# This is what is done to get the VARNA diagrams from Biers
				total_bpp = max(bpp_matrix[ii][cur_end], total_bpp)
				# print("%d %d %d %f\n" % (cur_cnt, ii, cur_end, bpp_matrix[ii][cur_end]))
				ii += 1
				cur_end -= 1
				continue

			if secstruct[ii] != '(':
				ii += 1
			elif cur_end != bps[ii]:
				cur_end -= 1

			if cur_cnt > 0:
				bpps += [total_bpp]
				cnts += [cur_cnt]

			cur_cnt = 0
			total_bpp = 0
		
		bpps += [total_bpp]
		cnts += [cur_cnt]

		#print(cnts)
		#print(bpps)
		# How many base-pairs are in stems passing the base-pair probability threshold?
		best_bpp = 0
		for ii, cnt in enumerate(cnts):
			#print(bpps[ii])
			if bpps[ii] > bpp_thresh:
				if cnt > len_cutoff: 
					best_bpp = max(bpps[ii], best_bpp)

		if best_bpp > 0:
			return True, best_bpp

		return False, -1
	
	def has_stem_dG(self, bp, seq, mfe, bpp_matrix=None, min_num_bp=6, bpp_cutoff=0.7):
		dG_lim = self.DG_LIM
		stems = self.get_stems(mfe)
		stems_region = []
		if not self.DOFIRST:
			stems_region = self.get_stems_region(stems, self.region, bp, len(seq)) # Second step
		else:
			stems_region = self.get_stems_region(stems, self.region, 0, bp) # First step

		has_zipper_stem = False
		best_dG = 200 # Some large number
		best_stem = ""
		for stem in stems_region:
			if stem[4] < min_num_bp:
				continue

			passes_bpp = True
			if bpp_matrix is not None:
				passes_bpp, best_bpp = self.get_best_bpp(mfe, stem, bpp_matrix, \
					bpp_thresh=bpp_cutoff)
			if not passes_bpp:
				continue

			dG = 200
			cur_stem = ""
			if not self.DOFIRST: 
				[dG, cur_stem] = self.get_max_stem_dG(stem, seq, bp + self.region[0], len(seq)-self.region[1]) # Second step
			else:
				[dG, cur_stem] = self.get_max_stem_dG(stem, seq, self.region[0], bp-self.region[1]) # First step
			if (dG < best_dG):
				has_zipper_stem = True
				best_dG = dG
				best_stem = cur_stem
		return [has_zipper_stem, best_stem, best_dG]

	def get_score_mfe(self, intron):
		[has_dG, best_stem, best_dG] = self.has_stem_dG(intron.bp, intron.seq, intron.mfe)
		if (has_dG):
			return best_dG
		else:
			return self.DG_LIM

	def get_score_mfe_bpp(self, intron):
		[has_dG, best_stem, best_dG] = self.has_stem_dG(intron.bp, intron.seq, \
			intron.mfe, bpp_matrix=intron.bpp)
		if (has_dG):
			return best_dG
		else:
			return self.DG_LIM

	def get_score_ens(self, intron):
		dg_total = 0
		dg_cnt = 0
		secstructs = intron.ens[:self.max_ens]
		for secstruct in secstructs:
			[has_dG, best_stem, best_dG] = self.has_stem_dG(intron.bp, intron.seq, secstruct)
			if has_dG:
				dg_total += best_dG
				dg_cnt += 1
		if dg_cnt > 0: 
			return dg_total/dg_cnt
		else:
			return 0
