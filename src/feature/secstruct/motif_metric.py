from feature.secstruct.secstruct_metric import SecstructMetric
import re
import math


class MotifMetric(SecstructMetric):
	def __init__(self, max_matches=100, name="MotifMetric"):
		self.max_matches = max_matches
		self.name = name

	# Expand characters into regex searchable versions of the character
	# 
	# char1 and char2 are treated as if they should be base paired if both
	# are passed in. 
	def expand_chars(self, char1, char2):
		# N is . for the regex pattern matching later.
		map_vars = {'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'], \
			'N': ['.'], 'R': ['A', 'G'], 'Y': ['C', 'T']}
		valid_bps = [('C', 'G'), ('G', 'C'), ('A', 'T'), ('T', 'A')]
		cgs = [('C', 'G'), ('G', 'C')]
		if char1 != "" and char2 != "":
			if char1 == 'M' and char2 == 'M':
				return valid_bps
			elif char1 == 'N' and char2 == 'N':
				return [('.', '.')]
			elif char1 == 'O' and char2 == 'O':
				return cgs
			else:
				return [(char1, char2)]
		elif char1 != "":
			return [(x, "") for x in map_vars[char1]]
		else:
			return [("", x) for x in map_vars[char2]]

	def add_variants(self, all_variants, char1, char2, pos1=-1, pos2=-1):
		new_variants = []
		options = self.expand_chars(char1, char2)
		for variant in all_variants:
			for option in options:
				new_variant = ""
				if (pos1 > -1):
					new_variant = variant[0:pos1] + option[0]
					if (pos2 > -1):
						new_variant += variant[pos1:-pos2] + option[1] + variant[-pos2:]
					else:
						new_variant += variant[pos1:] + option[1]
				else:
					new_variant = option[0]
					if (pos2 > -1):
						new_variant += variant[:-pos2] + option[1] + variant[-pos2:]
					else:
						new_variant += variant[:] + option[1]					
				new_variants += [new_variant]
		return new_variants

	def update_seqs(self, seq, bp, do_start):
		if do_start:
			return [seq[1:], bp[1:]]
		return [seq[:(len(seq)-1)], bp[:(len(bp)-1)]]

	# Create all possible regex patterns for this pair of sequences that matches the
	# supplied nucleotide sequence pair. 
	def make_all_patterns(self, seq, bp_constraint):
		all_variants = [""]
		pos1 = 0 # How many from the start to insert ?
		pos2 = 0 # How many from the end to insert ?
		while len(seq) > 0:
			char1 = seq[0]
			bp1 = bp_constraint[0]
			char2 = ""
			if (bp1 == '('):
				char2 = seq[len(seq)-1]
				bp2  = bp_constraint[len(bp_constraint)-1]
				while len(seq) > 0 and (bp2 != ')'):
					all_variants = self.add_variants(all_variants, "", char2, pos2=pos2)
					[seq, bp_constraint] = self.update_seqs(seq, bp_constraint, False)
					pos2 += 1
					if len(seq) > 0:
						char2 = seq[len(seq)-1]
						bp2 = bp_constraint[len(bp_constraint)-1]
				all_variants = self.add_variants(all_variants, char1, char2, pos1=pos1, pos2=pos2)
				[seq, bp_constraint] = self.update_seqs(seq, bp_constraint, True)
				[seq, bp_constraint] = self.update_seqs(seq, bp_constraint, False)
				pos1 += 1
				pos2 += 1
			else:
				all_variants = self.add_variants(all_variants, char1, "", pos1=pos1)
				[seq, bp_constraint] = self.update_seqs(seq, bp_constraint, True)
				pos1 += 1

		return all_variants

	def flip_constraint(self, constraint):
		constraint = constraint[::-1]
		constraint = constraint.replace('(', 'P')
		constraint = constraint.replace(')', '(')
		constraint = constraint.replace('P', ')')	
		return constraint


class MotifSingleMetric(MotifMetric):
	def __init__(self, seq, bp_cons, name="MotifSingleMetric"):
		super().__init__()
		self.seq = seq
		self.bp_cons = bp_cons
		self.all_patterns = self.make_all_patterns(seq, bp_cons)
		self.bp_cons_flip = self.flip_constraint(bp_cons)
		self.name = name

	def get_match_locs(self, intron_seq, do_rev=False):
		match_locs = []
		for pattern in self.all_patterns:
			if do_rev:
				pattern = pattern[::-1]
			for m in re.finditer(pattern, intron_seq):
				if (len(match_locs) < self.max_matches):
					match_locs += [m.start()]
		return match_locs
	
	def get_total_matches(self, match_locs, rev_match_locs, secstruct):
		total_matches = 0
		cons_len = len(self.bp_cons)
		for match in match_locs:
			if (secstruct[match:(match+cons_len)] == self.bp_cons):
				total_matches += 1
		for rev_match in rev_match_locs:
			if (secstruct[rev_match:(rev_match+cons_len)] == self.bp_cons_flip):
				total_matches += 1
		return total_matches

	def get_avg_total_matches(self, intron_seq, secstructs):
		match_locs = self.get_match_locs(intron_seq)
		rev_match_locs = self.get_match_locs(intron_seq, True)
		total_matches = 0
		for secstruct in secstructs:
			num_matches = self.get_total_matches(match_locs, rev_match_locs, secstruct)
			total_matches += num_matches
		return total_matches/len(secstructs)

	def get_score_mfe(self, intron):
		total_matches = self.get_avg_total_matches(intron.seq, [intron.mfe])
		if total_matches < 0.0001:
			total_matches = 0.0001
		return -math.log(total_matches)

	def get_score_ens(self, intron):
		total_matches = self.get_avg_total_matches(intron.seq, intron.ens)
		if total_matches < 0.0001:
			total_matches = 0.0001
		return -math.log(total_matches)


class MotifPairMetric(MotifMetric):
	def __init__(self, seq1, seq2, bp_cons1, bp_cons2, name="MotifPairMetric"):
		super().__init__()
		self.seq1 = seq1
		self.seq2 = seq2
		self.bp_cons1 = bp_cons1
		self.bp_cons2 = bp_cons2
		self.all_patterns = self.make_all_pair_patterns()
		self.bp_cons1_flip = self.flip_constraint(bp_cons1)
		self.bp_cons2_flip = self.flip_constraint(bp_cons2)
		self.name = name

	def make_all_pair_patterns(self):
		all_patterns = self.make_all_patterns(self.seq1 + self.seq2, self.bp_cons1 + self.bp_cons2)
		first_len = len(self.seq1)
		all_patterns = [(x[0:first_len], x[first_len:]) for x in all_patterns]
		return all_patterns

	# Get all potential location matches
	def get_match_locs(self, intron_seq, do_rev=False):
		match_pairs = []
		for pattern in self.all_patterns:
			pat0 = pattern[0]
			pat1 = pattern[1]
			if do_rev:
				tmp = pat0
				pat0 = pat1[::-1]
				pat1 = tmp[::-1] 
			match1_locs = []
			match2_locs = []
			for m in re.finditer(pat0, intron_seq):
				match1_locs += [m.start()]
			for m in re.finditer(pat1, intron_seq):
				match2_locs += [m.start()]
			for match1 in match1_locs:
				for match2 in match2_locs:
					match_dir_check = match1 + len(pat0) < match2
					if match_dir_check and len(match_pairs) < self.max_matches:
						match_pairs += [(match1, match2)]
		return match_pairs

	def get_total_matches(self, match_locs, rev_match_locs, secstruct):
		total_matches = 0
		len1 = len(self.bp_cons1)
		len2 = len(self.bp_cons2)
		for match in match_locs:
			if (secstruct[match[0]:(match[0]+len1)] == self.bp_cons1) and \
				(secstruct[match[1]:(match[1]+len2)] == self.bp_cons2):
				total_matches += 1
		for rev_match in rev_match_locs:
			if (secstruct[rev_match[0]:(rev_match[0]+len2)] == self.bp_cons2_flip) and \
				(secstruct[rev_match[1]:(rev_match[1]+len1)] == self.bp_cons1_flip):
				total_matches += 1
		return total_matches

	def get_avg_total_matches(self, intron_seq, secstructs):
		match_locs = self.get_match_locs(intron_seq)
		rev_match_locs = self.get_match_locs(intron_seq, True)
		total_matches = 0
		for secstruct in secstructs:
			num_matches = self.get_total_matches(match_locs, rev_match_locs, secstruct)
			total_matches += num_matches
		return total_matches/len(secstructs)

	def get_score_mfe(self, intron):
		return -self.get_avg_total_matches(intron.seq, [intron.mfe])

	def get_score_ens(self, intron):
		return -self.get_avg_total_matches(intron.seq, intron.ens)


class KinkturnMetric(MotifPairMetric):
	def __init__(self):
		super().__init__('MAG', 'AGNNRM', '(((', '))xxx)', name="KinkturnMetric")

class GNRAMetric(MotifSingleMetric):
	def __init__(self):
		super().__init__('MGNRAM', '(....)', name="GNRAMetric")
