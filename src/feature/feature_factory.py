from feature.secstruct.localization_metric import LocalizationMetric
from feature.secstruct.protection_metric import StartProtectionMetric, EndProtectionMetric, BPProtectionMetric
from feature.secstruct.motif_metric import KinkturnMetric, GNRAMetric
from feature.secstruct.stem_metric import StartToBPStemMetric, BPToEndStemMetric
from feature.secstruct.zipper_stem_metric import ZipperStemMetric

def get_feature_from_name(feature_name):
	if feature_name == "LocalizationMetric":
		return LocalizationMetric(pair1=0, pair2=-1, start_bp=False, end_bp=True)
	else:
		return eval(feature_name)()