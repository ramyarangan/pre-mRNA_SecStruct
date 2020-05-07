class Feature:
	def __init__(self):
		self.name = 'Feature'

	def get_full_name(self, feature_options):
		return self.name

	def apply(self, intron, feature_options):
		raise NotImplementedError