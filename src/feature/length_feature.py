from feature.feature import Feature

class LengthFeature(Feature):
	def __init__(self, name="LengthFeature"):
		self.name = name

	def apply(self, intron, feature_options):
		return len(intron.seq)