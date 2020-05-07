from util import features

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'mfe', 
					'verbose': True
					}

feature_df = features.get_features(['LocalizationMetric', 'StartToBPStemMetric'], 'test', 
	feature_options_all={'LocalizationMetric': secstruct_options, 'StartToBPStemMetric': secstruct_options})

print(feature_df.columns)