from feature.feature import Feature

class LengthFeature(Feature):
	def __init__(self, name="LengthFeature"):
		self.name = name

	def apply(self, intron, feature_options):
		start = intron.fivess_offset
		end = len(intron.seq) - intron.threess_offset
		return len(intron.seq[start:end])