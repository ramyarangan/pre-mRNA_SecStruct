from util import features_db

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'mfe', 
					'verbose': True,
					'force_eval': False
					}

all_features = ["ZipperStemMetric", "LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", \
		"StartProtectionMetric", "EndProtectionMetric", "BPProtectionMetric"]#, \
		# "HasEarlyStopFeature", "ThreeprimeDistStopFeature", "RPKMFeature"]


feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

feature_df = features_db.get_features(all_features, 'standard_shuffle', feature_options_all=feature_options_all)
print(feature_df.columns)