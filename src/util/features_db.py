from core.intron import IntronSet
from feature.feature_factory import get_feature_from_name
from util.sec_struct import add_secstruct_mfe_to_database, add_secstruct_ens_to_database
from feature.secstruct.secstruct_metric import SecstructMetric
from config import DATABASE_PATH
import os
import pandas as pd

def get_features_filename(intron_class):
	return DATABASE_PATH + 'introns/' + intron_class + '/features.csv'

def build_intron_set(intron_class, feature_options={}):
	intron_seq_file = DATABASE_PATH + 'introns/' + intron_class + '/base_info.dat' 
	if not os.path.isfile(intron_seq_file):
		raise RuntimeError('Intron class missing base_info.dat file')

	mfe_file = ""
	ens_file = ""
	if 'secstruct_pkg' in feature_options.keys():
		sec_struct_pref = feature_options['secstruct_pkg']
		mfe_file_check = DATABASE_PATH + 'introns/' + intron_class + '/' + sec_struct_pref + '_mfe.dat'
		ens_file_check = DATABASE_PATH + 'introns/' + intron_class + '/' + sec_struct_pref + '_ens.dat'
		if os.path.isfile(mfe_file_check):
			mfe_file = mfe_file_check
		if os.path.isfile(ens_file_check):
			ens_file = ens_file_check
	
	all_introns = IntronSet()
	all_introns.init_from_files(intron_seq_file, mfe_filename=mfe_file, ens_filename=ens_file)

	return all_introns

def setup_features_file_from_base_info(intron_class, all_introns):
	features_file = get_features_filename(intron_class)

	col_names = ['ChrPos', 'Strand', 'Name']
	df = pd.DataFrame(columns=col_names)

	for ii in range(len(all_introns.introns)):
		intron = all_introns.introns[ii]
		new_entry = [intron.name]
		df.loc[len(df)] = [intron.chr_pos, intron.strand, intron.name]

	write_file = open(features_file, 'w')
	df.to_csv(path_or_buf=write_file, index=False)
	write_file.close()

def get_feature_full_name(feature_name, feature_options):
	feature = get_feature_from_name(feature_name)
	return feature.get_full_name(feature_options)

def get_features_df(intron_class, feature_names=[], feature_options_all={}):
	features_file = get_features_filename(intron_class)
	df_file = open(features_file)
	features_df = pd.read_csv(df_file)
	df_file.close()

	if feature_names != []:
		full_names = [get_feature_full_name(feature_name, feature_options_all[feature_name]) for \
			feature_name in feature_names]
		features_df = features_df[full_names]

	return features_df

# Requires that all introns have non-na features for this value, otherwise reevaluate
# Can change later.
def get_features_in_database(intron_class):
	features_df = get_features_df(intron_class)
	columns = features_df.columns
	no_na_columns = []
	for column in columns:
		if not features_df[column].isnull().values.any():
			no_na_columns += [column]
	return no_na_columns

def get_feature_vals(feature, all_introns, feature_options):
	verbose = False
	if 'verbose' in feature_options.keys():
		verbose = feature_options['verbose']

	feature_vals = []
	for ii, intron in enumerate(all_introns.introns):
		if verbose:
			print("Feature: %s, Evaluating intron: %d of %d" % \
				(feature.name, ii, len(all_introns.introns)))
		feature_vals += [feature.apply(intron, feature_options)]
		# feature_vals += [1] # Fast evaluation for now

	return feature_vals

# Needs to also return all_introns because it can get updated to include
# the secondary structure here.
def get_feature_vals_update_secstruct(feature_name, intron_class, feature_options, all_introns):
	feature = get_feature_from_name(feature_name)
	
	# Make sure that intron set has secondary structures computed if needed
	if isinstance(feature, SecstructMetric):
		if 'secstruct_pkg' not in feature_options.keys():
			raise RuntimeError("Feature options must include the secondary structure package info")
		if 'secstruct_type' not in feature_options.keys():
			raise RuntimeError("Feature options must include the secondary structure estimation type (mfe or ens)")

		sec_struct_pref = feature_options['secstruct_pkg']
		sec_struct_type = feature_options['secstruct_type']

		struct_file = DATABASE_PATH + 'introns/' + intron_class + '/' + sec_struct_pref + '_' + sec_struct_type + '.dat'
		if not os.path.isfile(struct_file):
			if sec_struct_type == 'mfe':
				add_secstruct_mfe_to_database(intron_class, sec_struct_pref)
			if sec_struct_type == 'ens':
				add_secstruct_ens_to_database(intron_class, sec_struct_pref)

		all_introns = build_intron_set(intron_class, feature_options)

	feature_vals = get_feature_vals(feature, all_introns, feature_options)
	return [feature_vals, all_introns]

def add_features_to_database(features_to_add, intron_class, feature_options_all, all_introns):
	features_df = get_features_df(intron_class)
	
	for feature_name in features_to_add:
		feature_options = feature_options_all[feature_name]
		full_feature_name = get_feature_full_name(feature_name, feature_options)
		[feature_vals, all_introns] = get_feature_vals_update_secstruct(feature_name, intron_class, \
			feature_options, all_introns)
		features_df[full_feature_name] = feature_vals
	
	features_file = get_features_filename(intron_class)	
	write_file = open(features_file, 'w')
	features_df.to_csv(path_or_buf=write_file, index=False)
	write_file.close()

# FeatureData holds options for features: secondary structure package, ens/mfe, print progress or no. 
# It is a dictionary from feature name to feature options
def get_features(feature_names, intron_class, feature_options_all={}):
	all_introns = build_intron_set(intron_class)

	features_file = get_features_filename(intron_class)
	if not os.path.isfile(features_file):
		setup_features_file_from_base_info(intron_class, all_introns)
	
	# Get the names of features that are stored in the database
	features_in_database = get_features_in_database(intron_class)

	# Add any features that are requested that are missing in the database
	# Important: this might change all_introns as more data is loaded when needed.
	features_to_add = []
	for feature_name in feature_names:
		feature_options = feature_options_all[feature_name]
		feature_full_name = get_feature_full_name(feature_name, feature_options)
		if feature_full_name not in features_in_database:
			features_to_add += [feature_name]
		elif 'force_eval' in feature_options and feature_options['force_eval']:
			features_to_add += [feature_name]
	if len(features_to_add) > 0:
		add_features_to_database(features_to_add, intron_class, feature_options_all, all_introns)

	# Get a dataframe of the features requested
	return get_features_df(intron_class, feature_names=feature_names, feature_options_all=feature_options_all)