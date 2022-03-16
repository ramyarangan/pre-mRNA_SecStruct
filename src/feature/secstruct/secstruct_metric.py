from feature.feature import Feature

class SecstructMetric(Feature):
	def __init__(self):
		self.name = 'SecstructMetric'

	def get_full_name(self, feature_options):
		if 'secstruct_pkg' not in feature_options:
			raise RuntimeError("SecstructMetric feature requires secondary structure package in feature_options")
		if 'secstruct_type' not in feature_options:
			raise RuntimeError("SectructMetric feature requires secondary structure type: ens or mfe")
		full_name = self.name + '_' + feature_options['secstruct_pkg'] + '_' + feature_options['secstruct_type']
		if 'use_bpp' in feature_options and feature_options['use_bpp']:
			full_name += '_' + "BPP"
		return full_name

	def apply(self, intron, feature_options):
		if 'secstruct_type' not in feature_options:
			raise RuntimeError("SectructMetric feature requires secondary structure type: ens or mfe")
		
		if feature_options['secstruct_type'] == 'mfe':
			if 'use_bpp' in feature_options and feature_options['use_bpp']:
				return get_score_mfe_bpp(intron)
			return self.get_score_mfe(intron)

		if feature_options['secstruct_type'] == 'ens':
			if 'use_bpp' in feature_options and feature_options['use_bpp']:
				raise RuntimeError("Cannot use BPP matrix when making ensemble predictions")
			return self.get_score_ens(intron)
	
	def get_score_mfe_bpp(self, intron):
		return self.get_score_mfe(intron)

	def get_score_mfe(self, intron):
		raise NotImplementedError

	def get_score_ens(self, intron):
		raise NotImplementedError
