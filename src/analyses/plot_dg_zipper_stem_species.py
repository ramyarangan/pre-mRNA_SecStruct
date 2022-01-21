"""
Make plots comparing intron and control zipper stem free energies across Saccharomyces species
Example usage: 
 python analyses/plot_dg_zipper_stem_species.py standard_min_50_max_600 standard_min_50_max_600_shuffle --make_violin
"""
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd 
import numpy as np 
import argparse
import os 

from util import features_db
from config import DATABASE_PATH

parser = argparse.ArgumentParser(description='Parameters for plotting dG stem features across species')
parser.add_argument('intron_class', type=str, help='Intron class')
parser.add_argument('control_class', type=str, help='Control class')
parser.add_argument('--make_histograms', default=False, action='store_true', \
	help='Make histograms for intron and control histograms of dG values for each species; print counts and bins')
parser.add_argument('--make_violin', default=False, action='store_true', \
	help='Make violin plots compraing intron and control dG values across all species')
args = parser.parse_args()

intron_class = args.intron_class
control_class = args.control_class
make_histograms = args.make_histograms
make_violin = args.make_violin
intron_class = 'species_hooks/' + intron_class
control_class = 'species_hooks/' + control_class 

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'ens', 
					'verbose': True,
					'force_eval': False
					}

all_features = ["ZipperStemStartMetric", "ZipperStemEndMetric"]

metric_names = [features_db.get_feature_full_name(feature_name, secstruct_options) \
		for feature_name in all_features]

feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

intron_path = os.path.join(DATABASE_PATH, 'introns/' + intron_class)

species_names = []
for species_name in os.listdir(intron_path):
	if species_name == 'scer':
		continue
	species_names += [species_name]

df_dict = {}
for feature in all_features:
	df_dict[feature] = pd.DataFrame(columns=['species', 'dG_diff'])

for species_name in species_names:
	print("Species: %s\n" % species_name)
	if species_name == 'scer':
		continue
	intron_class_species = intron_class + '/' + species_name
	standard_feature_df = features_db.get_features(all_features, \
		intron_class_species, feature_options_all=feature_options_all)
	standard_feature_df = standard_feature_df.dropna(axis=0)

	control_class_species = control_class + '/' + species_name
	control_feature_df = features_db.get_features(all_features, \
		control_class_species, feature_options_all=feature_options_all)
	control_feature_df = control_feature_df.dropna(axis=0)

	if make_histograms:
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

	for ii, metric in enumerate(metric_names):
		dg_diff_vals = []
		for idx, row in standard_feature_df.iterrows():
			dg_diff_vals += [row[metric]]
		for idx, row in control_feature_df.iterrows():
			dg_diff_vals[idx] -= row[metric]

		cur_df = df_dict[all_features[ii]] 
		for dg_diff_val in dg_diff_vals:
			cur_df.loc[len(cur_df)] = [species_name, dg_diff_val]

if make_violin:
	for feature in all_features:
		plt.figure(figsize=(8,2))
		ax = sns.violinplot(data=df_dict[feature], x="species", y="dG_diff", cut=0, scale='width')
		ax.axhline(0, color='black', linestyle='--')
		plt.ylabel("Intron - Control dG of zipper stems")
		plt.xticks(fontsize=8)
		plt.yticks(np.arange(-25, 11, 10), fontsize=8)
		plt.show()

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
