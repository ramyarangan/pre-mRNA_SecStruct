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
import seaborn as sns
import pandas as pd

secstruct_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'mfe', 
					'verbose': True,
					'force_eval': False
					}

all_features = ["ZipperStemMetric", "LocalizationMetric", "StartToBPStemMetric", "BPToEndStemMetric", \
		"StartProtectionMetric", "EndProtectionMetric"]

all_features = ["ZipperStemMetric",  "HasEarlyStopFeature", "ThreeprimeDistStopFeature",]

all_features = ["ZipperStemMetric"]# "ThreeprimeDistStopFeature", "RPKMFeature"]


feature_options_all = {}
for feature in all_features:
	feature_options_all[feature] = secstruct_options

standard_feature_df = features_db.get_features(all_features, 'standard', feature_options_all=feature_options_all)
standard_feature_df = standard_feature_df.dropna(axis=0)
decoy_feature_df = features_db.get_features(all_features, 'decoy', feature_options_all=feature_options_all)
decoy_feature_df = decoy_feature_df.dropna(axis=0)

all_feature_df = pd.concat([standard_feature_df, decoy_feature_df], axis=0)
exp_class_data = np.array([0] * standard_feature_df.shape[0] + [1] * decoy_feature_df.shape[0])
# np.random.shuffle(exp_class_data)
feature_matrix = all_feature_df.to_numpy()
print(all_feature_df.isna().sum())

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
	for cutoff in bins:
		pred_pos = np.zeros(len(pred_probs))
		pred_pos[pred_probs > cutoff] = 1
		tpr = sum(pred_pos * exp_classes)/sum(exp_classes)
		fpr = sum(pred_pos * (1 - exp_classes))/sum(1 - exp_classes)
		tprs += [tpr]
		fprs += [fpr]
	
	auc = sklearn.metrics.auc(fprs, tprs)
	print("AUC: %f" % auc)
	
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

		# Train decision tree 
		clf = ensemble.RandomForestClassifier(n_estimators=n_estimators,
										 max_leaf_nodes=max_leaf_nodes, 
										 min_samples_leaf=min_sample_leaf)
		clf = clf.fit(feature_matrix[train_idxs,:], exp_data[train_idxs])

		# Predict on test data
		if len(test_idxs) > 0:
			predictions = clf.predict_proba(feature_matrix[test_idxs,:])
			all_pred_data += list(predictions)
			all_exp_data += list(exp_data[test_idxs])

		# Predict on train data
		if len(train_idxs) > 0:
			train_predictions = clf.predict_proba(feature_matrix[train_idxs,:])
			all_pred_data_train += list(train_predictions)
			all_exp_data_train += list(exp_data[train_idxs])
	return (all_pred_data, all_exp_data, all_pred_data_train, all_exp_data_train)

# Monte-Carlo Cross Validation

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