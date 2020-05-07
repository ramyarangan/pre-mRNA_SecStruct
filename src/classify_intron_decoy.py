from util import features

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'mfe', 
					'verbose': True
					}

#all_features = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", "StartProtectionMetric", 
#		"EndProtectionMetric", "BPProtectionMetric"]
all_features = ["ZipperStemMetric"]

feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

feature_df = features.get_features(all_features, 'standard', feature_options_all=feature_options_all)
print(feature_df.columns)