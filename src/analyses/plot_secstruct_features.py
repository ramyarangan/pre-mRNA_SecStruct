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

from util import features_db
from config import DATABASE_PATH

PLOT_NAME_DICT = {
	"LocalizationMetric": "5'SS to Branchpoint Distance",
	"StartToBPStemMetric": "5'SS to Branchpoint Stem",
	"BPToEndStemMetric": "Branchpoint to 3'SS Stem", 
	"StartProtectionMetric": "5'SS Protection", 
	"EndProtectionMetric": "3'SS Protection", 
	"BPProtectionMetric": "Branchpoint Protection"
}

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
	for species_name in os.listdir(intron_path):
		all_species += [species_name]

	df = pd.DataFrame(index=np.array(plot_names), columns=np.array(all_species))
	
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

		for ii, metric in enumerate(metric_names):
			intron_vals = []
			control_vals = []
			for idx, row in standard_feature_df.iterrows():
				intron_vals += [row[metric]]
			for idx, row in control_feature_df.iterrows():
				control_vals += [row[metric]]
			_, pval = stats.wilcoxon(np.array(intron_vals), np.array(control_vals))
			df.at[plot_names[ii], species_name] = pval
	df = np.log(df.astype(float))

	plt.figure(figsize=(10,8))
	# sns.heatmap(df, cmap='coolwarm', center=0.3, annot=True, vmin=0, vmax=0.5)
	sns.heatmap(df, cmap='coolwarm', center=0.5, annot=False, vmin=-10, vmax=0)
	plt.xticks(rotation=45, fontsize=8)
	plt.show()

def make_violin_plots(standard_feature_df, \
	control_feature_df, metric_names, plot_names):
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

	max_vals = {}
	min_vals = {}
	for metric in metric_names: 
		max_intron = max(standard_feature_df[metric])
		max_control = max(control_feature_df[metric])
		max_vals[metric] = max(max_intron, max_control)

		min_intron = min(standard_feature_df[metric])
		min_control = min(control_feature_df[metric])
		min_vals[metric] = min(min_intron, min_control)

	df = pd.DataFrame(columns=plot_names)
	for idx, row in standard_feature_df.iterrows():
		new_entry = []
		for metric in metric_names:
			norm_intron = (row[metric]-min_vals[metric])/(max_vals[metric]-min_vals[metric])
			control_row = control_feature_df.iloc[[idx]]
			norm_control = (control_row[metric]-min_vals[metric])/(max_vals[metric]-min_vals[metric])
			new_entry += [norm_intron - norm_control]
		df.loc[len(df)] = new_entry

	plt.figure(figsize=(10,1.5))
	ax = sns.violinplot(data=df)# , inner="stick")
	ax.axhline(0, color='black', linestyle='--')
	plt.ylabel("Intron Metric - Control Metric")
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
						'verbose': False,
						'force_eval': False
						}

	all_features = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", "StartProtectionMetric", 
			"EndProtectionMetric", "BPProtectionMetric"]

	metric_names = [features_db.get_feature_full_name(feature_name, secstruct_options) \
			for feature_name in all_features]

	plot_names = [PLOT_NAME_DICT[x] for x in all_features]

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

		make_violin_plots(standard_feature_df, control_feature_df, metric_names, plot_names)

	if make_species_heatmap:
		make_heatmap(all_features, feature_options_all, \
			intron_class_species, control_class_species, plot_names)

