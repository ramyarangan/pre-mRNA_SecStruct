import pandas as pd 
import sys
from scipy import stats 
import numpy as np
from matplotlib import pyplot as plt
from util import features_db

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'ens', 
					'verbose': False,
					'force_eval': False
					}

all_features = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", "StartProtectionMetric", 
		"EndProtectionMetric", "BPProtectionMetric"]

metric_names = [features_db.get_feature_full_name(feature_name, secstruct_options) \
		for feature_name in all_features]

feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

standard_feature_df = features_db.get_features(all_features, 'standard_extend20', feature_options_all=feature_options_all)
standard_feature_df = standard_feature_df.dropna(axis=0)
shift_feature_df = features_db.get_features(all_features, 'standard_extend20_shift', feature_options_all=feature_options_all)
shift_feature_df = shift_feature_df.dropna(axis=0)

ii = 1
for metric in metric_names:
	intron_vals = []
	control_vals = []
	for idx, row in standard_feature_df.iterrows():
		intron_vals += [row[metric]]
	for idx, row in shift_feature_df.iterrows():
		control_vals += [row[metric]]
	plt.subplot(1, len(metric_names), ii)
	ii += 1
	plt.hist(np.array(intron_vals) - np.array(control_vals), alpha=0.5, rwidth=0.85, color="forestgreen")
	plt.axvline(x=0, linestyle='--', color='black')
	print(metric)
	print(stats.ttest_rel(np.array(intron_vals), np.array(control_vals)))
	print(stats.wilcoxon(np.array(intron_vals), np.array(control_vals)))

plt.show()
