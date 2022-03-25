from feature.secstruct.secstruct_metric import SecstructMetric

class ProtectionMetric(SecstructMetric):
	def __init__(self, length, do_bp=False, do_end=False, name="ProtectionMetric"):
		self.length = length
		self.do_bp = do_bp
		self.do_end = do_end
		self.name = name

	def get_avg_protection(self, secstructs, start_idx, bp, intron_len):
		if (self.do_bp):
			start_idx = bp - 4
		if (self.do_end):
			start_idx = intron_len - self.length
		total_protection = 0
		for secstruct in secstructs:
			cur_protection = 0
			max_idx = min(start_idx + self.length, len(secstruct))
			for ii in range(start_idx, max_idx):
				if (secstruct[ii] != '.'):
					cur_protection += 1
			cur_protection = cur_protection/(max_idx - start_idx)
			total_protection += cur_protection
		return total_protection/len(secstructs)

	def get_score_mfe(self, intron):
		return self.get_avg_protection([intron.mfe], intron.fivess_offset, 
			intron.bp, len(intron.seq) - intron.threess_offset)

	def get_score_ens(self, intron):
		return self.get_avg_protection(intron.ens, intron.fivess_offset, 
			intron.bp, len(intron.seq) - intron.threess_offset)


class StartProtectionMetric(ProtectionMetric):
	def __init__(self):
		super().__init__(7, name="StartProtectionMetric")


class EndProtectionMetric(ProtectionMetric):
	def __init__(self):
		super().__init__(7, do_end=True, name="EndProtectionMetric")


class BPProtectionMetric(ProtectionMetric):
	def __init__(self):
		super().__init__(7, do_bp=True, name="BPProtectionMetric")