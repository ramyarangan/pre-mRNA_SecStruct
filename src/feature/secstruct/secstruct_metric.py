from feature.feature import Feature

class SecstructMetric(Feature):
	def __init__(self):
		self.name = 'SecstructMetric'

	def get_full_name(self, feature_options):
		if 'secstruct_pkg' not in feature_options:
			raise RuntimeError("SecstructMetric feature requires secondary structure package in feature_options")
		if 'secstruct_type' not in feature_options:
			raise RuntimeError("SectructMetric feature requires secondary structure type: ens or mfe")
		return self.name + '_' + feature_options['secstruct_pkg'] + '_' + feature_options['secstruct_type']

	def apply(self, intron, feature_options):
		if 'secstruct_type' not in feature_options:
			raise RuntimeError("SectructMetric feature requires secondary structure type: ens or mfe")
		
		if feature_options['secstruct_type'] == 'mfe':
			return self.get_score_mfe(intron)

		if feature_options['secstruct_type'] == 'ens':
			return self.get_score_ens(intron)
	
	def get_score_mfe(self, intron):
		raise NotImplementedError

	def get_score_ens(self, intron):
		raise NotImplementedError
