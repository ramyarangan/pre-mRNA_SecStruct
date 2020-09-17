from feature.secstruct.localization_metric import LocalizationMetric
from feature.secstruct.protection_metric import StartProtectionMetric, EndProtectionMetric, BPProtectionMetric
from feature.secstruct.motif_metric import KinkturnMetric, GNRAMetric
from feature.secstruct.stem_metric import StartToBPStemMetric, BPToEndStemMetric
from feature.secstruct.zipper_stem_metric import ZipperStemMetric
from feature.early_stop_feature import HasEarlyStopFeature, ThreeprimeDistStopFeature
from feature.rpkm_feature import RPKMFeature

def get_feature_from_name(feature_name):
	if feature_name == "LocalizationMetric":
		return LocalizationMetric(start_bp=False, end_bp=True)
	else:
		return eval(feature_name)()