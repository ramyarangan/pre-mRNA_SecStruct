from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd 
import numpy as np 
import sys

from util import features_db


intron_class = sys.argv[1] # E.g. standard_allsize_min_50_max_600
control_class = sys.argv[2] # E.g. standard_allsize_min_50_max_600_shuffle

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'ens', 
					'verbose': True,
					'force_eval': False
					}

all_features = ["ZipperStemStartMetric", "ZipperStemEndMetric"]

metric_names = [features_db.get_feature_full_name(feature_name, secstruct_options) \
		for feature_name in all_features]

# Determined by looking at the data distribution
plot_mins = {
	metric_names[0]: -25, 
	metric_names[1]: -17.5
}

feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

standard_feature_df = features_db.get_features(all_features, \
	intron_class, feature_options_all=feature_options_all)
standard_feature_df = standard_feature_df.dropna(axis=0)
control_feature_df = features_db.get_features(all_features, \
	control_class, feature_options_all=feature_options_all)
control_feature_df = control_feature_df.dropna(axis=0)

for ii, feature in enumerate(all_features):
	mask = standard_feature_df[metric_names[ii]] < -0.05
	counts, bins, bars = plt.hist(standard_feature_df[mask][metric_names[ii]], \
		range = [plot_mins[metric_names[ii]], 0])
	plt.show()
	print(counts)
	print(bins)

	mask = control_feature_df[metric_names[ii]] < -0.05
	counts, bins, bars = plt.hist(control_feature_df[mask][metric_names[ii]], \
		range = [plot_mins[metric_names[ii]], 0])
	plt.show()
	print(counts)
	print(bins)

"""
[ 1.  0.  2.  9. 18. 29. 36. 38. 33. 19.]
[-25.  -22.5 -20.  -17.5 -15.  -12.5 -10.   -7.5  -5.   -2.5   0. ]
[ 0.  0.  0.  0.  2.  5.  7. 24. 42. 46.]
[-25.  -22.5 -20.  -17.5 -15.  -12.5 -10.   -7.5  -5.   -2.5   0. ]
[ 1.  0.  0.  1.  4.  5.  6. 14. 25. 33.]
[-25.  -22.5 -20.  -17.5 -15.  -12.5 -10.   -7.5  -5.   -2.5   0. ]
[ 0.  0.  0.  0.  0.  1.  9.  9. 16. 49.]
[-25.  -22.5 -20.  -17.5 -15.  -12.5 -10.   -7.5  -5.   -2.5   0. ]
"""
