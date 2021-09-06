import pandas as pd 
import sys
from scipy import stats 
import numpy as np
from matplotlib import pyplot as plt
from util import features_db
import sys
import os

from config import DATABASE_PATH

intron_class = sys.argv[1] # E.g. standard_allsize_min_50_max_600
control_class = sys.argv[2] # E.g. standard_allsize_min_50_max_600_shuffle

intron_class = 'species_hooks/' + intron_class
control_class = 'species_hooks/' + control_class 

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'ens', 
					'verbose': True,
					'force_eval': False
					}

	
all_features = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", "StartProtectionMetric", 
		"EndProtectionMetric", "BPProtectionMetric"]

metric_names = [features_db.get_feature_full_name(feature_name, secstruct_options) \
		for feature_name in all_features]

feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

intron_path = os.path.join(DATABASE_PATH, 'introns/' + intron_class)
for species_name in os.listdir(intron_path):
	print("Species: %s\n" % species_name)

	intron_class_species = intron_class + '/' + species_name
	standard_feature_df = features_db.get_features(all_features, \
		intron_class_species, feature_options_all=feature_options_all)
	standard_feature_df = standard_feature_df.dropna(axis=0)

	control_class_species = control_class + '/' + species_name
	control_feature_df = features_db.get_features(all_features, \
		control_class_species, feature_options_all=feature_options_all)
	control_feature_df = control_feature_df.dropna(axis=0)

	# standard_feature_df = features_db.get_features(all_features, \
	#  	'standard', feature_options_all=feature_options_all)
	# standard_feature_df = standard_feature_df.dropna(axis=0)
	# control_feature_df = features_db.get_features(all_features, \
	# 	'standard_shuffle', feature_options_all=feature_options_all)
	# control_feature_df = control_feature_df.dropna(axis=0)

	print(standard_feature_df.shape)
	print(control_feature_df.shape)

	ii = 1
	for metric in metric_names:
		intron_vals = []
		control_vals = []
		for idx, row in standard_feature_df.iterrows():
			intron_vals += [row[metric]]
		for idx, row in control_feature_df.iterrows():
			control_vals += [row[metric]]
		#plt.subplot(1, len(metric_names), ii)
		ii += 1
		#plt.hist(np.array(intron_vals) - np.array(control_vals), alpha=0.5, rwidth=0.85, color="forestgreen")
		#plt.axvline(x=0, linestyle='--', color='black')
		print(metric)
		print(stats.ttest_rel(np.array(intron_vals), np.array(control_vals)))
		print(stats.wilcoxon(np.array(intron_vals), np.array(control_vals)))

	#plt.show()
