from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd 
import numpy as np 
import sys

from util import features_db
from util.features_db import * 
from util.gene_names import *

intron_class = sys.argv[1] # E.g. standard_allsize_min_50_max_600
control_class = sys.argv[2] # E.g. standard_allsize_min_50_max_600_shuffle

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'ens', 
					'verbose': True,
					'force_eval': False
					}

all_features = ["ZipperStemStartMetric", "ZipperStemEndMetric"]

summary_cutoffs = {
	"ZipperStemStartMetric": -10, 
	"ZipperStemEndMetric": -10
}

metric_names = [features_db.get_feature_full_name(feature_name, secstruct_options) \
		for feature_name in all_features]

# Determined by looking at the data distribution
plot_mins = {
	metric_names[0]: -45, # -25, 
	metric_names[1]: -45# -17.5
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

# Print number of introns 
check_update_secstruct(intron_class, secstruct_options)
all_introns = build_intron_set(intron_class, secstruct_options)

rpg_mask = get_rpg_intron_mask(all_introns)
print("Total number of RPGs: %d" % sum(rpg_mask))

intron_len_mask = [(len(x.seq) > 200) for x in all_introns.introns]
intron_len_mask = np.array(intron_len_mask)
print("Total introns with length > 200: %d" % sum(intron_len_mask))

if len(standard_feature_df) != len(all_introns.introns):
	print("Don't have the expected number of features")

# Get number of stems that pass threshold 
for ii, feature in enumerate(all_features):
	total_standard = len(standard_feature_df)
	total_control = len(control_feature_df)

	summary_cutoff_val = summary_cutoffs[feature]
	
	print(min(standard_feature_df[metric_names[ii]]))
	mask = standard_feature_df[metric_names[ii]] < summary_cutoff_val
	mask_len = np.logical_and(np.array(mask), intron_len_mask)
	num_passing_standard = sum(mask)
	num_passing_standard_len = sum(mask_len)
	num_passing_standard_rpg = sum(np.logical_and(np.array(mask), rpg_mask))

	mask = control_feature_df[metric_names[ii]] < summary_cutoff_val
	mask_len = np.logical_and(np.array(mask), intron_len_mask)
	num_passing_control = sum(mask)
	num_passing_control_len = sum(mask_len)
	num_passing_control_rpg = sum(np.logical_and(np.array(mask), rpg_mask))

	print(feature)
	print("Num passing cutoff %f kcal/mol - introns: %d/%d" % \
		(summary_cutoff_val, num_passing_standard, total_standard))
	print("Num passing cutoff %f kcal/mol - controls: %d/%d" % \
		(summary_cutoff_val, num_passing_control, total_control))
	print("Num passing cutoff %f kcal/mol with len > 200 - introns: %d/%d" % \
		(summary_cutoff_val, num_passing_standard_len, sum(intron_len_mask)))
	print("Num passing cutoff %f kcal/mol with len > 200 - controls: %d/%d" % \
		(summary_cutoff_val, num_passing_control_len, sum(intron_len_mask)))
	print("Num passing cutoff %f kcal/mol in RPGs - introns: %d/%d" % \
		(summary_cutoff_val, num_passing_standard_rpg, sum(rpg_mask)))
	print("Num passing cutoff %f kcal/mol in RPGs - controls: %d/%d" % \
		(summary_cutoff_val, num_passing_control_rpg, sum(rpg_mask)))

# Get histogram bin values - used for plotting in GraphPad Prism
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
