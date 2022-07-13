"""
Make plots showing pvalues from secondary structure feature comparisons 

Either plot violin plots for all features for a specific intron/control pair, or make
a p-value heatmap for all species

Example usage: 
python analyses/plot_secstruct_features.py --make_species_heatmap --intron_class_species standard_min_50_max_600 --control_class_species standard_min_50_max_600_shuffle_seq_matched
"""
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd 
import numpy as np 
import argparse
import os
from scipy import stats 
import math
from util import features_db
from config import DATABASE_PATH

PLOT_NAME_DICT = {
	"LocalizationMetric": "5'SS-BP Dist",
	"StartToBPStemMetric": "5'SS-BP Stem",
	"BPToEndStemMetric": "BP-3'SS Stem", 
	"StartProtectionMetric": "5'SS Protection", 
	"EndProtectionMetric": "3'SS Protection", 
	"BPProtectionMetric": "BP Protection", 
	"ZipperStemStartMetric": "Zipper dG",
	"ZipperStemEndMetric": "End dG",
	"LongestStemMetric": "Longest Stem", 
	"NWJMetric": "# NWJs", 
	"MLDMetric": "MLD",
	"LengthFeature": "Length"
}
SPECIES_NAMES = ['scer', 'smik', 'skud', 'suva', 'cgla', 'kafr', 'knag', 'ncas', \
	'ndai', 'tbla', 'tpha', 'kpol', 'zrou', 'tdel', 'klac', 'agos', 'ecym', 'sklu', \
	'kthe', 'kwal']

def read_pvals_from_file(filename):
	species_dir = {}
	f = open(filename)
	lines = f.readlines()
	f.close()

	cur_species = ""
	count = 0
	for line in lines:
		if line == "\n":
			continue
		if ((line.split()[0] != "Comparing") and (line.split()[0] != "Intron")) and \
			(line.split()[0] != "p-val") and (line.split()[0] != "No"):
			cur_species = line.split()[0]
		if (line.split()[0] == "p-val"):
			cur_line = line[0:(len(line)-2)]
			count += 1
			if count % 2 == 0:
				new_pval = float(cur_line.split()[5])
				new_pval = min(new_pval, 1 - new_pval)
				if cur_species in species_dir:
					species_dir[cur_species] += [new_pval]
				else:
					species_dir[cur_species] = [new_pval]

	return species_dir

def make_heatmap_from_species_dir(species_dir):
	metrics = ["Localization", "5'BPStem", "BP3'Stem", "StartProtection", "EndProtection", "BPProtection"]
	ordered_names = ["scer", "smik", "skud", "suva", "cgla", "kafr", \
		"knag", "ncas", "ndai", "tbla", "tpha", "kpol", "zrou", "tdel", "klac", \
		"agos", "ecym", "sklu", "kthe", "kwal"]
	array = []
	for name in ordered_names:
		array += [species_dir[name]]
		print(name)
		print(len(species_dir[name]))
	array = np.array(array)
	array = array[:,[0, 2, 3, 4, 5, 6]]
	print(array)

	plt.figure(figsize=(10,8))
	sns.heatmap(df, cmap='coolwarm', center=0.3, annot=True, vmin=0, vmax=0.5)
	# sns.heatmap(df, cmap='coolwarm', center=0.5, annot=False, vmin=0, vmax=1)
	plt.show()

def make_heatmap(all_features, feature_options_all, \
	intron_class, control_class, plot_names):
	intron_class = 'species_hooks/' + intron_class
	control_class = 'species_hooks/' + control_class 
	intron_path = os.path.join(DATABASE_PATH, 'introns/' + intron_class)

	all_species = []
	for species_name in SPECIES_NAMES:
		if species_name == 'scer':
			continue
		all_species += [species_name]

	df = pd.DataFrame(index=np.array(plot_names), columns=np.array(all_species))
	
	for species_name in SPECIES_NAMES:
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

		for ii, metric in enumerate(metric_names):
			intron_vals = []
			control_vals = []
			for idx, row in standard_feature_df.iterrows():
				intron_vals += [row[metric]]
			for idx, row in control_feature_df.iterrows():
				control_vals += [row[metric]]
			_, pval = stats.wilcoxon(np.array(intron_vals), np.array(control_vals))
			df.at[plot_names[ii], species_name] = pval
	df = np.log10(df.astype(float))

	plt.figure(figsize=(10,8))
	# sns.heatmap(df, cmap='coolwarm', center=0.3, annot=True, vmin=0, vmax=0.5)
	sns.heatmap(df, cmap='coolwarm', center=0.5, annot=False, vmin=-4, vmax=0)
	plt.xticks(rotation=45, fontsize=8)
	plt.show()

def get_percentile(arr, perc):
	arr = np.array(arr)
	return np.percentile(arr, perc * 100)

def make_violin(standard_feature_df, \
	control_feature_df, metric_names, plot_names, all_introns, do_ratio_norm):
	print(standard_feature_df.shape)
	print(control_feature_df.shape)

	# norm_long_stem = []
	# norm_long_stem_control = []
	# longest_stem_vals = standard_feature_df["LongestStemMetric_Vienna_ens"]
	# longest_stem_vals_controls = control_feature_df["LongestStemMetric_Vienna_ens"]
	# for ii, intron in enumerate(all_introns.introns):
	# 	norm_long_stem += [longest_stem_vals[ii]/len(intron.seq)]
	# 	norm_long_stem_control += [longest_stem_vals_controls[ii]/len(intron.seq)]
	# standard_feature_df["LongestStemMetric_Vienna_ens"] = norm_long_stem
	# control_feature_df["LongestStemMetric_Vienna_ens"] = norm_long_stem_control

	max_vals = {}
	min_vals = {}
	perc_95_vals = {}
	perc_5_vals = {}
	for metric in metric_names: 
		all_vals = list(standard_feature_df[metric])
		all_vals += list(control_feature_df[metric])

		# if "LongestStem" in metric:
		# 	plt.hist(list(standard_feature_df[metric]), alpha=0.5, label='intron')
		# 	plt.hist(list(control_feature_df[metric]), alpha=0.5, label='control')
		# 	plt.legend(loc='upper right')
		# 	plt.show()
		max_vals[metric] = max(all_vals)
		perc_95_vals[metric] =  get_percentile(all_vals, 0.95)

		min_vals[metric] = min(all_vals)
		perc_5_vals[metric] =  get_percentile(all_vals, 0.05)
		print("%s: %f-%f\n" % (metric, perc_5_vals[metric], perc_95_vals[metric]))

		intron_vals = standard_feature_df[metric]
		control_vals = control_feature_df[metric]
		
		print(stats.ttest_rel(np.array(intron_vals), np.array(control_vals)))
		print(stats.wilcoxon(np.array(intron_vals), np.array(control_vals)))

	df = pd.DataFrame(columns=plot_names)
	for idx, row in standard_feature_df.iterrows():
		new_entry = []
		for ii, metric in enumerate(metric_names):
			do_ratio = do_ratio_norm[ii]
			perc_range = perc_95_vals[metric]-perc_5_vals[metric]
			min_max_range = max_vals[metric] - min_vals[metric]
			norm_intron = (row[metric]-min_vals[metric])/min_max_range
			control_row = control_feature_df.iloc[[idx]]
			norm_control = (control_row[metric]-min_vals[metric])/min_max_range
			norm_diff = norm_intron - norm_control
			norm_ratio = 2 * norm_diff/(norm_control + norm_intron)
			if do_ratio:
				new_entry += [norm_ratio]
			else:
				new_entry += [norm_diff]
		df.loc[len(df)] = new_entry

	plt.figure(figsize=(4,8))
	ax = sns.violinplot(data=df, scale="count")# , inner="stick")
	ax.axhline(0, color='black', linestyle='--')
	plt.ylabel("Intron Metric - Control Metric")
	plt.ylim(-3.5, 3.5)
	plt.xticks(rotation=45)
	plt.show()

# def make_countplot(df_metric_filename):
# 	f = open(df_metric_filename)
# 	features = pd.read_csv(f)
# 	f.close()
# 
# 	df = pd.DataFrame(columns=["Category","More Efficient"])
# 
# 	metrics = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", "StartProtectionMetric", 
# 		"EndProtectionMetric", "BPProtectionMetric"] # ,  "GNRAMetric"
# 	new_names = ["5'SS to Branchpoint Distance", "5'SS to Branchpoint Stem", "Branchpoint to 3'SS Stem", \
# 		"5'SS Protection", "3'SS Protection", "Branchpoint Protection"] # , "GNRA Motif Presence"
# 	rename_cols = {}
# 	for ii in range(len(metrics)):
# 		rename_cols[metrics[ii]] = new_names[ii]
# 
# 	for idx, row in features.iterrows():
# 		for metric in metrics: 
# 			if (row[metric] < row["Control_"+metric]):
# 				df.loc[len(df)] = [rename_cols[metric], "Intron"]
# 			elif (row[metric] > row["Control_"+metric]):
# 				df.loc[len(df)] = [rename_cols[metric], "Control"]
# 			else:
# 				df.loc[len(df)] = [rename_cols[metric], "Tie"]
# 
# 
# 
# 	plt.figure(figsize=(20,3))
# 	sns.set_palette("muted")
# 	ax = sns.countplot(data=df, x="Category", hue="More Efficient", hue_order=["Intron", "Control", "Tie"],\
# 		palette={"Intron": "C0", "Control": "C1", "Tie": "C2"})
# 	ax.legend_.remove()
# 	plt.ylabel("Count")
# 	plt.show()

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Parameters for plotting secondary structure feature figures')
	parser.add_argument('--make_violin_plots', default=False, action='store_true', \
		 help='Make violin plots comparing an intron set and control set')
	parser.add_argument('--intron_class_violin', type=str, \
		help='Intron class with features computed; e.g. standard_allsize_min_50_max_600')
	parser.add_argument('--control_class_violin', type=str, \
		help='Control class with features computed; e.g. standard_allsize_min_50_max_600_shuffle')
	parser.add_argument('--make_species_heatmap', default=False, action='store_true', \
		 help='Make violin plots comparing an intron set and control set')
	parser.add_argument('--intron_class_species', type=str, \
		help='Directory to species intron classes with features computed; e.g. species_hooks/standard_min_50_max_600')
	parser.add_argument('--control_class_species', type=str, \
		help='Directory to species control classes with features computed; e.g. species_hooks/standard_min_50_max_600_phylo')
	args = parser.parse_args()

	make_violin_plots = args.make_violin_plots
	intron_class_violin = ""
	control_class_violin = ""
	if make_violin_plots:
		intron_class_violin = args.intron_class_violin
		control_class_violin = args.control_class_violin
	make_species_heatmap = args.make_species_heatmap
	if make_species_heatmap:
		intron_class_species = args.intron_class_species
		control_class_species = args.control_class_species

	secstruct_options = {'secstruct_pkg': 'Vienna', 
						'secstruct_type': 'ens', 
						'use_bpp': False,
						'verbose': False,
						'force_eval': False
						}


	all_features = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", "StartProtectionMetric", 
		"EndProtectionMetric", "BPProtectionMetric", "ZipperStemStartMetric", "ZipperStemEndMetric", \
		"LongestStemMetric", "NWJMetric", "MLDMetric"]
	all_features = ["ZipperStemStartMetric", "ZipperStemEndMetric", \
		"LongestStemMetric", "MLDMetric", "LocalizationMetric", "LengthFeature"]
	do_ratio_norm = [True, True, True, True, True]

	metric_names = [features_db.get_feature_full_name(feature_name, secstruct_options) \
			for feature_name in all_features]
	all_metric_names = metric_names
	metric_names = metric_names[:-1]

	plot_names = [PLOT_NAME_DICT[x] for x in all_features]
	plot_names = plot_names[:-1]

	feature_options_all = {}
	for feature in all_features:
		feature_options_all[feature] = secstruct_options

	if make_violin_plots:
		standard_feature_df = features_db.get_features(all_features, \
			intron_class_violin, feature_options_all=feature_options_all)
		standard_feature_df = standard_feature_df.dropna(axis=0)
		control_feature_df = features_db.get_features(all_features, \
			control_class_violin, feature_options_all=feature_options_all)
		control_feature_df = control_feature_df.dropna(axis=0)

		all_introns = features_db.build_intron_set(intron_class_violin)
		make_violin(standard_feature_df, control_feature_df, metric_names, plot_names, all_introns, do_ratio_norm)

	if make_species_heatmap:
		make_heatmap(all_features, feature_options_all, \
			intron_class_species, control_class_species, plot_names)

