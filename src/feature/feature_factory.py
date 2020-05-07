from secstruct.localization_metric import LocalizationMetric
from secstruct.protection_metric import StartProtectionMetric, EndProtectionMetric, BPProtectionMetric
from secstruct.motif_metric import KinkturnMetric, GNRAMetric
from secstruct.stem_metric import StartToBPStemMetric, BPToEndStemMetric

def get_feature_from_name(feature_name)
	if feature_name == "LocalizationMetric":
		return LocalizationMetric(pair1=0, pair2=-1, start_bp=False, end_bp=True)
	else
		return eval(feature_name)()