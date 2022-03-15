from feature.secstruct.localization_metric import LocalizationMetric
from feature.secstruct.protection_metric import StartProtectionMetric, EndProtectionMetric, BPProtectionMetric
from feature.secstruct.motif_metric import KinkturnMetric, GNRAMetric
from feature.secstruct.stem_metric import StartToBPStemMetric, BPToEndStemMetric
from feature.secstruct.zipper_stem_metric import ZipperStemMetric
from feature.secstruct.longest_stem_metric import LongestStemMetric
from feature.secstruct.nwj_metric import NWJMetric
from feature.secstruct.mld_metric import MLDMetric
from feature.early_stop_feature import HasEarlyStopFeature, ThreeprimeDistStopFeature
from feature.rpkm_feature import RPKMFeature

def get_feature_from_name(feature_name):
	if feature_name == "LocalizationMetric":
		return LocalizationMetric(start_bp=False, end_bp=True)
	if feature_name == "ZipperStemStartMetric":
		return ZipperStemMetric(do_first=True, name="ZipperStemStartMetric")
	if feature_name == "ZipperStemEndMetric":
		return ZipperStemMetric(do_first=False, name="ZipperStemEndMetric")
	else:
		return eval(feature_name)()