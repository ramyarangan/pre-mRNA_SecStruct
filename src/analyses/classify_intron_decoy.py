from util import features_db
import numpy as np
from sklearn import linear_model
from sklearn import tree
from sklearn import ensemble
from sklearn import metrics
import sklearn
import shap
import argparse
import os
from matplotlib import pyplot as plt
from scipy.stats import spearmanr
from imblearn.over_sampling import RandomOverSampler
import seaborn as sns
import pandas as pd
import random
from collections import Counter

SECSTRUCT_FEATURES = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", \
			"StartProtectionMetric", "EndProtectionMetric", "BPProtectionMetric", "ZipperStemStartMetric", \
			"ZipperStemEndMetric", "LongestStemMetric", "MLDMetric"]

# SECSTRUCT_FEATURES = ["LocalizationMetric", "StartProtectionMetric", "EndProtectionMetric",
#	"BPProtectionMetric", "ZipperStemStartMetric", "ZipperStemEndMetric", "LongestStemMetric", "MLDMetric"]

SECSTRUCT_FEATURES = ["LocalizationMetric", "ZipperStemStartMetric", "ZipperStemEndMetric", "LongestStemMetric", "MLDMetric"]

# SECSTRUCT_FEATURES = ["ZipperStemStartMetric", "ZipperStemEndMetric", "LongestStemMetric", "MLDMetric"]

def filter_decoys(df, columns_list, cutoff_list, cutoff_is_greater, do_control=False, control_df=None):
	for ii, column in enumerate(columns_list):
		cutoff = cutoff_list[ii]
		idxs = (df[column] <= cutoff)
		if cutoff_is_greater[ii]:
			idxs = (df[column] > cutoff)
		if do_control:
			control_df = control_df[idxs]
		df = df[idxs]
	return df, control_df

def update_feature_list(df, feature_list):
	feature_list = [features_db.get_feature_full_name(feature_name, secstruct_options) \
		for feature_name in feature_list]
	df = df[feature_list]
	return df

def expore_cutoffs():
	cutoff_val = 2
	cutoff_val = 50
	print(np.sum(standard_feature_df.to_numpy()[:,0] > cutoff_val)/standard_feature_df.shape[0])
	print(np.sum(standard_feature_df.to_numpy()[:,0] < cutoff_val)/standard_feature_df.shape[0])
	print(np.sum(decoy_feature_df.to_numpy()[:,0] > cutoff_val)/decoy_feature_df.shape[0])
	print(np.sum(decoy_feature_df.to_numpy()[:,0] < cutoff_val)/decoy_feature_df.shape[0])

def print_tree(feature_matrix, exp_class_data, feature_names):
	clf = tree.DecisionTreeClassifier(max_leaf_nodes = 5, min_samples_leaf = 3)
	clf = clf.fit(feature_matrix, exp_class_data)
	r = tree.export_text(clf, feature_names=feature_names)
	print(r)
	# ThreeprimeDistStopFeature <= 2 for standard introns
	# RPKMFeature >= 582.74 for standard introns

def get_auc(pred_probs, exp_classes, bin_width=0.01):
	pred_probs = np.array(pred_probs)
	exp_classes = np.array(exp_classes)
	bins = np.arange(0, 1 + bin_width, bin_width)
	tprs = []
	fprs = []
	for cutoff in bins:
		pred_pos = np.zeros(len(pred_probs))
		pred_pos[pred_probs > cutoff] = 1
		tpr = sum(pred_pos * exp_classes)/sum(exp_classes)
		fpr = sum(pred_pos * (1 - exp_classes))/sum(1 - exp_classes)
		tprs += [tpr]
		fprs += [fpr]
	
	auc = sklearn.metrics.auc(fprs, tprs)
	return auc

def get_fprs_tprs(pred_probs, exp_classes, bin_width):
	pred_probs = np.array(pred_probs)
	exp_classes = np.array(exp_classes)
	bins = np.arange(0, 1 + bin_width, bin_width)
	tprs = []
	fprs = []
	cohen_kappas = []
	for cutoff in bins:
		pred_pos = np.zeros(len(pred_probs))
		pred_pos[pred_probs > cutoff] = 1
		cohen_kappa = sklearn.metrics.cohen_kappa_score(pred_pos, exp_classes)
		cohen_kappas += [cohen_kappa]
		tpr = sum(pred_pos * exp_classes)/sum(exp_classes)
		fpr = sum(pred_pos * (1 - exp_classes))/sum(1 - exp_classes)
		tprs += [tpr]
		fprs += [fpr]
	return [fprs, tprs, cohen_kappas, bins]

def calc_auc(pred_probs, exp_classes, bin_width=0.01):
	[fprs, tprs, _, _] = \
		get_fprs_tprs(pred_probs, exp_classes, bin_width)
	auc = sklearn.metrics.auc(fprs, tprs)
	print("AUC: %f" % auc)
	return auc

def plot_auc_curve(pred_probs, exp_classes, bin_width=0.01):
	[fprs, tprs, cohen_kappas, bins] = \
		get_fprs_tprs(pred_probs, exp_classes, bin_width)
	
	plt.plot(bins, cohen_kappas)
	plt.xlabel("False Postive Rate")
	plt.ylabel("Cohen's Kappa")
	plt.show()

	plt.plot(fprs, tprs, color='black')
	plt.plot(bins, bins, color='red', linestyle='--')
	plt.xlabel("False Positive Rate")
	plt.ylabel("True Positive Rate")
	plt.title("ROC curve for intron/ decoy classification")
	plt.show()

def get_train_test_idxs(num_items, num_train):
	rnd_idxs = np.random.permutation(num_items)
	train_idxs = rnd_idxs[0:num_train]
	test_idxs = rnd_idxs[(num_train+1):]
	return (train_idxs, test_idxs)

def random_forest_clf_mccv(feature_matrix, exp_data, 
					   num_train, niter, n_estimators=10, max_leaf_nodes=10, min_sample_leaf=5, 
					   do_shap=False):
	all_pred_data = []
	all_exp_data = []

	all_pred_data_train = []
	all_exp_data_train = []

	exp_data = np.array(exp_data)
	
	# Train many models
	for ii in range(niter):
		print("%d of %d iterations" % (ii, niter))

		# Randomize order of training data
		(train_idxs, test_idxs) = get_train_test_idxs(np.size(feature_matrix, 0), num_train)


		ros = RandomOverSampler(random_state=42)
		# print(Counter(exp_data[train_idxs]))
		feature_train, data_train = ros.fit_resample(feature_matrix[train_idxs,:], exp_data[train_idxs])
		# print(Counter(data_train))

		# Train decision tree 
		clf = ensemble.RandomForestClassifier(n_estimators=n_estimators,
										 max_leaf_nodes=max_leaf_nodes, 
										 min_samples_leaf=min_sample_leaf)
		clf = clf.fit(feature_train, data_train)

		# Predict on test data
		if len(test_idxs) > 0:
			predictions = clf.predict_proba(feature_matrix[test_idxs,:])
			all_pred_data += list(predictions)
			all_exp_data += list(exp_data[test_idxs])

		# Predict on train data
		if len(train_idxs) > 0:
			train_predictions = clf.predict_proba(feature_train)
			all_pred_data_train += list(train_predictions)
			all_exp_data_train += list(data_train)

		if do_shap:
			explainer = shap.TreeExplainer(clf)
			shap_values = explainer.shap_values(feature_train)
			plot_names = [x.replace("Metric", "") for x in SECSTRUCT_FEATURES]
			shap.summary_plot(shap_values[0], feature_train, feature_names=plot_names)
			plt.show()

	return (all_pred_data, all_exp_data, all_pred_data_train, all_exp_data_train)

def print_trees(feature_matrix, exp_data, column_names, niter, train_size=0.9):
	num_train = int(train_size * np.size(feature_matrix, 0))
	for ii in range(niter):
		(train_idxs, test_idxs) = get_train_test_idxs(np.size(feature_matrix, 0), num_train)
		ros = RandomOverSampler(random_state=42)
		feature_train, data_train = ros.fit_resample(feature_matrix[train_idxs,:], exp_data[train_idxs])
		print_tree(feature_train, data_train, column_names)

# Monte-Carlo Cross Validation
def do_mc_cv_random_forest(feature_matrix, exp_class_data, do_plots=False, do_shap=False):
	# Training data size
	train_size = 0.9
	num_train = int(train_size * np.size(feature_matrix, 0))
	niter = 10 # 100
	if do_shap:	
		niter = 10
	n_estimators = 5000
	max_leaf_nodes = 5
	min_sample_leaf = 3

	[all_pred_data, all_exp_data, all_pred_data_train, all_exp_data_train] = \
		random_forest_clf_mccv(feature_matrix, exp_class_data, num_train, niter,
						   max_leaf_nodes=max_leaf_nodes, min_sample_leaf=min_sample_leaf, 
						  n_estimators=n_estimators, do_shap=do_shap)

	print("Training data")
	pred_probs = np.array([x[1] for x in all_pred_data_train])
	calc_auc(pred_probs, all_exp_data_train)
	if do_plots:
		plot_auc_curve(pred_probs, all_exp_data_train)

	print("Test data")
	pred_probs = np.array([x[1] for x in all_pred_data])
	calc_auc(pred_probs, all_exp_data)
	if do_plots:
		plot_auc_curve(pred_probs, all_exp_data)

def get_class_control_diff_df(standard_feature_df, control_feature_df):
	max_vals = {}
	min_vals = {}
	for metric in metric_names: 
		all_vals = list(standard_feature_df[metric])
		all_vals += list(control_feature_df[metric])
		max_vals[metric] = max(all_vals)
		min_vals[metric] = min(all_vals)

	df = pd.DataFrame(columns=standard_feature_df.columns)
	for idx, row in standard_feature_df.iterrows():
		new_entry = []
		for metric in metric_names:
			min_max_range = max_vals[metric] - min_vals[metric]
			norm_intron = (row[metric]-min_vals[metric])/min_max_range
			control_row = control_feature_df.iloc[[idx]]
			norm_control = (control_row[metric]-min_vals[metric])/min_max_range
			new_entry += [norm_intron - norm_control]
		df.loc[len(df)] = new_entry
	return df

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Parameters for classifying intron sets')
	parser.add_argument('intron_class', type=str, help='Intron class to classify from decoys')
	parser.add_argument('decoy_class', type=str, help='Decoy class to classify from introns')
	parser.add_argument('--intron_control_class', default='', type=str, help='Intron control class')
	parser.add_argument('--decoy_control_class', default='', type=str, help='Decoy control class')
	parser.add_argument('--do_enrichment', default=False, action='store_true', \
		 help='Classify by enrichment')
	parser.add_argument('--do_dropout', default=False, action='store_true', \
		 help='Do dropout on each feature')
	parser.add_argument('--do_print_trees', default=False, action='store_true', \
		 help='Print some sample decision trees')
	parser.add_argument('--do_roc_plot', default=False, action='store_true', \
		 help='Plot Cohen kappa curves and AUC curves for train and test set')
	parser.add_argument('--do_shap', default=False, action='store_true', \
		 help='Do Shap feature importance analysis and plot Shap scores')
	args = parser.parse_args()

	intron_class = args.intron_class
	decoy_class = args.decoy_class
	intron_control_class = args.intron_control_class
	decoy_control_class = args.decoy_control_class
	do_enrichment = args.do_enrichment
	do_dropout = args.do_dropout
	do_print_trees = args.do_print_trees
	do_plots = args.do_roc_plot
	do_shap = args.do_shap

	secstruct_options = {'secstruct_pkg': 'RNAstructure_DMS', 
					'secstruct_type': 'mfe', 
					'use_bpp': True,
					'verbose': False,
					'force_eval': False
					}

	secstruct_features = SECSTRUCT_FEATURES

	all_features = secstruct_features + ["LengthFeature"] # ["ThreeprimeDistStopFeature", "RPKMFeature", "LengthFeature"] # "HasEarlyStopFeature"

	feature_options_all = {}
	for feature in all_features:
		feature_options_all[feature] = secstruct_options

	standard_feature_df = features_db.get_features(all_features, intron_class, feature_options_all=feature_options_all)
	standard_feature_df = standard_feature_df.dropna(axis=0)
	decoy_feature_df = features_db.get_features(all_features, decoy_class, feature_options_all=feature_options_all)
	decoy_feature_df = decoy_feature_df.dropna(axis=0)
	standard_control_feature_df = None
	decoy_control_feature_df = None
	if do_enrichment:
		standard_control_feature_df = features_db.get_features(all_features, \
			intron_control_class, feature_options_all=feature_options_all)
		standard_control_feature_df = standard_control_feature_df.dropna(axis=0)
		decoy_control_feature_df = features_db.get_features(all_features, \
			decoy_control_class, feature_options_all=feature_options_all)
		decoy_control_feature_df = decoy_control_feature_df.dropna(axis=0)

	decoy_feature_df, decoy_control_feature_df = filter_decoys(decoy_feature_df, 
		["LengthFeature"], [150], [True], do_control=do_enrichment, control_df=decoy_control_feature_df)
	#decoy_feature_df, decoy_control_feature_df = filter_decoys(decoy_feature_df, 
	#	["ThreeprimeDistStopFeature", "RPKMFeature", "LengthFeature"], 
	#	[50, 50, 150], [False, True, True], do_control=do_enrichment, control_df=decoy_control_feature_df)
	standard_feature_df, standard_control_feature_df = filter_decoys(standard_feature_df, 
		["LengthFeature"], [150], [True], do_control=do_enrichment, control_df=standard_control_feature_df)
	print(standard_feature_df.shape)
	print(decoy_feature_df.shape)

	feature_list = secstruct_features
	standard_feature_df = update_feature_list(standard_feature_df, feature_list)
	decoy_feature_df = update_feature_list(decoy_feature_df, feature_list)
	if do_enrichment:
		standard_control_feature_df = update_feature_list(standard_control_feature_df, feature_list)
		decoy_control_feature_df = update_feature_list(decoy_control_feature_df, feature_list)

	#plt.hist(standard_feature_df['LongestStemMetric_Vienna_ens'], alpha=0.5, label='intron')
	#plt.hist(decoy_feature_df['LongestStemMetric_Vienna_ens'], alpha=0.5, label='decoy')
	#plt.legend(loc='upper right')
	#plt.show()

	if do_enrichment:
		standard_feature_df = get_class_control_diff_df(standard_feature_df, standard_control_feature_df)
		decoy_feature_df = get_class_control_diff_df(decoy_feature_df, decoy_control_feature_df)

	all_feature_df = pd.concat([standard_feature_df, decoy_feature_df], axis=0)
	exp_class_data = np.array([0] * standard_feature_df.shape[0] + [1] * decoy_feature_df.shape[0])
	# np.random.shuffle(exp_class_data)
	feature_matrix = all_feature_df.to_numpy()

	do_mc_cv_random_forest(feature_matrix, exp_class_data, do_plots=do_plots, do_shap=do_shap)

	if do_print_trees:
		print_trees(feature_matrix, exp_class_data, list(all_feature_df.columns), 10)

	if do_dropout:
		for feature in secstruct_features:
			print("Dropout: %s\n" % feature)

			dropout_feature_list = []
			for tmp_feature in secstruct_features:
				if tmp_feature == feature:
					continue
				dropout_feature_list += [tmp_feature]
			dropout_feature_list = [feature]
			standard_dropout_df = \
				update_feature_list(standard_feature_df, dropout_feature_list)
			decoy_dropout_df = \
				update_feature_list(decoy_feature_df, dropout_feature_list)
			all_feature_df = pd.concat([standard_dropout_df, decoy_dropout_df], axis=0)
			exp_class_data = np.array([0] * standard_dropout_df.shape[0] + [1] * decoy_dropout_df.shape[0])
			# np.random.shuffle(exp_class_data)
			feature_matrix = all_feature_df.to_numpy()

			do_mc_cv_random_forest(feature_matrix, exp_class_data, do_plots=do_plots, do_shap=do_shap)
