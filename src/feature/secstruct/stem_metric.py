from feature.secstruct.secstruct_metric import SecstructMetric
import math 

class StemMetric(SecstructMetric):

	def __init__(self, do_bp_start=False, do_bp_end=False, \
		do_intron_end=False, dist_cutoff=-1, req_stem_len=-1, name="Stem Metric"):
		self.do_bp_start = do_bp_start
		self.do_bp_end = do_bp_end
		self.dist_cutoff = dist_cutoff
		self.req_stem_len = req_stem_len
		self.do_intron_end = do_intron_end
		self.name = name

	# Get base pairs
	def get_base_pairs(self, secstruct):
		bp_map = {}
		bp_stack = []
		for ii in range(0, len(secstruct)):
			if (secstruct[ii] == "("):
				bp_stack = bp_stack + [ii]
			if (secstruct[ii] == ")"):
				bp_start = bp_stack[-1]
				bp_stack = bp_stack[0:-1]
				bp_map[bp_start] = ii
				bp_map[ii] = bp_start
		return bp_map

	# maybe allow for bulges in this stem?
	def find_stem(self, bp_map, intron):
		start = intron.fivess_offset
		end = len(intron.seq) - intron.threess_offset
		if (self.do_bp_start):
			start = intron.bp
		if (self.do_bp_end):
			end = intron.bp
		if (self.do_intron_end):
			end = len(intron.seq) - intron.threess_offset
		end_thresh = max(end - self.dist_cutoff, start)
		start_thresh = min(end, start + self.dist_cutoff)
		stem_found = False
		for ii in range(start, start_thresh):
			if (ii in bp_map) and \
				(bp_map[ii] > end_thresh) and \
				(bp_map[ii] < end):
				cur_start = ii
				cur_end = bp_map[ii]
				num_bps = 1
				while ((cur_start + 1) in bp_map) and \
					(bp_map[cur_start + 1] == cur_end - 1):
					num_bps += 1
					cur_start += 1
					cur_end -= 1
				if (num_bps > self.req_stem_len):
					stem_found = True
			if (stem_found):
				break
		if stem_found:
			return 1
		return 0

	def get_score_mfe(self, intron):
		bp_map = self.get_base_pairs(intron.mfe)
		return -self.find_stem(bp_map, intron)

	def get_score_ens(self, intron):
		num_stems = 0
		for secstruct in intron.ens:
			bp_map = self.get_base_pairs(secstruct)
			num_stems += self.find_stem(bp_map, intron)
		avg_num_stems = num_stems/len(intron.ens)
		if avg_num_stems < 0.00001:
			avg_num_stems = 0.00001
		return -math.log(avg_num_stems)


class StartToBPStemMetric(StemMetric):
	def __init__(self):
		super().__init__(do_bp_end=True, dist_cutoff=15, req_stem_len=3, name="StartToBPStemMetric")


class BPToEndStemMetric(StemMetric):
	def __init__(self):
		super().__init__(do_bp_start=True, do_intron_end=True, dist_cutoff=15, req_stem_len=3, name="BPToEndStemMetric")

