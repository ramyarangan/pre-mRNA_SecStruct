from feature.feature import Feature
from util.general import get_rpkms
import numpy as np

class RPKMFeature(Feature):
	def __init__(self, name="RPKMFeature"):
		self.name = name
		self.rpkms = get_rpkms()

	def apply(self, intron, feature_options):
		if str(intron.ensembl_name) == 'nan':
			return np.nan
		if intron.ensembl_name not in self.rpkms:
			return np.nan
		return self.rpkms[intron.ensembl_name]