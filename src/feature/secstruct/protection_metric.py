from intron_metric import IntronMetric

class ProtectionMetric(IntronMetric):
	def __init__(self, length, start_idx=-1, do_bp=False, do_end=False, name="ProtectionMetric"):
		self.start_idx = start_idx
		self.length = length
		self.do_bp = do_bp
		self.do_end = do_end
		self.name = name

	def get_avg_protection(self, secstructs, bp, intron_len):
		start_idx = self.start_idx
		if (self.do_bp):
			start_idx = bp - 4
		if (self.do_end):
			start_idx = intron_len - 7
		total_protection = 0
		for secstruct in secstructs:
			cur_protection = 0
			for ii in range(start_idx, start_idx + self.length):
				if (secstruct[ii] != '.'):
					cur_protection += 1
			cur_protection = cur_protection/self.length
			total_protection += cur_protection
		return total_protection/len(secstructs)

	def get_score_mfe(self, intron):
		return self.get_avg_protection([intron.mfe], intron.bp, len(intron.seq))

	def get_score_ens(self, intron):
		return self.get_avg_protection(intron.ens, intron.bp, len(intron.seq))


class StartProtectionMetric(ProtectionMetric):
	def __init__(self):
		super().__init__(7, start_idx=0, name="StartProtectionMetric")


class EndProtectionMetric(ProtectionMetric):
	def __init__(self):
		super().__init__(7, do_end=True, name="EndProtectionMetric")


class BPProtectionMetric(ProtectionMetric):
	def __init__(self):
		super().__init__(7, do_bp=True, name="BPProtectionMetric")