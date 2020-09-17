from util import features_db
import numpy as np
from sklearn import linear_model
from sklearn import tree
from sklearn import ensemble
from sklearn import metrics
import sklearn
import os
from matplotlib import pyplot as plt
from scipy.stats import spearmanr
from imblearn.over_sampling import RandomOverSampler
import seaborn as sns
import pandas as pd
import random
from collections import Counter

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'ens', 
					'verbose': False,
					'force_eval': False
					}

secstruct_features = ["LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", \
		"StartProtectionMetric", "EndProtectionMetric"] # "ZipperStemMetric", 

all_features = secstruct_features + ["ThreeprimeDistStopFeature", "RPKMFeature"] # "HasEarlyStopFeature"

feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

def filter_decoys(df, columns_list, cutoff_list, cutoff_is_greater):
	for ii, column in enumerate(columns_list):
		cutoff = cutoff_list[ii]
		idxs = (df[column] <= cutoff)
		if cutoff_is_greater[ii]:
			idxs = (df[column] > cutoff)
		df = df[idxs]
	return df

def update_feature_list(df, feature_list):
	feature_list = [features_db.get_feature_full_name(feature_name, secstruct_options) \
		for feature_name in feature_list]
	df = df[feature_list]
	return df

standard_feature_df = features_db.get_features(all_features, 'standard', feature_options_all=feature_options_all)
standard_feature_df = standard_feature_df.dropna(axis=0)
decoy_feature_df = features_db.get_features(all_features, 'decoy', feature_options_all=feature_options_all)
decoy_feature_df = decoy_feature_df.dropna(axis=0)
# Filter out low transcription genes and NMD introns
decoy_feature_df = filter_decoys(decoy_feature_df, ["ThreeprimeDistStopFeature", "RPKMFeature"], [0, 50], [False, True])
print(decoy_feature_df.shape)

feature_list = secstruct_features
standard_feature_df = update_feature_list(standard_feature_df, feature_list)
decoy_feature_df = update_feature_list(decoy_feature_df, feature_list)

all_feature_df = pd.concat([standard_feature_df, decoy_feature_df], axis=0)
exp_class_data = np.array([0] * standard_feature_df.shape[0] + [1] * decoy_feature_df.shape[0])
# np.random.shuffle(exp_class_data)
feature_matrix = all_feature_df.to_numpy()



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

def plot_auc_curve(pred_probs, exp_classes, bin_width=0.01):
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
	
	auc = sklearn.metrics.auc(fprs, tprs)
	print("AUC: %f" % auc)
	
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
					   num_train, niter, n_estimators=10, max_leaf_nodes=10, min_sample_leaf=5):
	all_pred_data = []
	all_exp_data = []

	all_pred_data_train = []
	all_exp_data_train = []

	exp_data = np.array(exp_data)
	
	# Train many models
	for ii in range(niter):
		# Randomize order of training data
		(train_idxs, test_idxs) = get_train_test_idxs(np.size(feature_matrix, 0), num_train)


		ros = RandomOverSampler(random_state=42)
		print(Counter(exp_data[train_idxs]))
		feature_train, data_train = ros.fit_resample(feature_matrix[train_idxs,:], exp_data[train_idxs])
		print(Counter(data_train))

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
	return (all_pred_data, all_exp_data, all_pred_data_train, all_exp_data_train)

def print_trees(feature_matrix, exp_data, column_names, niter, train_size=0.9):
	num_train = int(train_size * np.size(feature_matrix, 0))
	for ii in range(niter):
		(train_idxs, test_idxs) = get_train_test_idxs(np.size(feature_matrix, 0), num_train)
		ros = RandomOverSampler(random_state=42)
		feature_train, data_train = ros.fit_resample(feature_matrix[train_idxs,:], exp_data[train_idxs])
		print_tree(feature_train, data_train, column_names)

# Monte-Carlo Cross Validation
def do_mc_cv_random_forest(feature_matrix, exp_class_data):
	# Training data size
	train_size = 0.9
	num_train = int(train_size * np.size(feature_matrix, 0))
	niter = 100
	n_estimators = 5
	max_leaf_nodes = 5
	min_sample_leaf = 3

	[all_pred_data, all_exp_data, all_pred_data_train, all_exp_data_train] = \
		random_forest_clf_mccv(feature_matrix, exp_class_data, num_train, niter,
						   max_leaf_nodes=max_leaf_nodes, min_sample_leaf=min_sample_leaf, 
						  n_estimators=n_estimators)

	print("Training data")
	pred_probs = np.array([x[1] for x in all_pred_data_train])
	plot_auc_curve(pred_probs, all_exp_data_train)

	print("Test data")
	pred_probs = np.array([x[1] for x in all_pred_data])
	plot_auc_curve(pred_probs, all_exp_data)

# do_mc_cv_random_forest(feature_matrix, exp_class_data)
print_trees(feature_matrix, exp_class_data, list(all_feature_df.columns), 10)

