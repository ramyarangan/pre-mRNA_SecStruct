from util import features_db

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'ens', 
					'verbose': True,
					'force_eval': True
					}

all_features = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", \
		"StartProtectionMetric", "EndProtectionMetric", "BPProtectionMetric", \
		"ZipperStemMetric", "HasEarlyStopFeature", "ThreeprimeDistStopFeature"]

feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

feature_df = features_db.get_features(all_features, 'decoy', feature_options_all=feature_options_all)
print(feature_df.columns)