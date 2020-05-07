from util import features_db

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'mfe', 
					'verbose': True,
					'force_eval': True
					}

#all_features = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", "StartProtectionMetric", 
#		"EndProtectionMetric", "BPProtectionMetric", "ZipperStemMetric", ]
all_features = ["HasEarlyStopFeature", "ThreeprimeDistStopFeature"]

feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

feature_df = features.get_features(all_features, 'standard', feature_options_all=feature_options_all)
print(feature_df.columns)