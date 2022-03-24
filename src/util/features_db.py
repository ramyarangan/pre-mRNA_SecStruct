import os
import pandas as pd

from core.intron import IntronSet
from feature.feature_factory import get_feature_from_name
from util.secstruct_db import add_secstruct_mfe_to_database, add_secstruct_ens_to_database
from feature.secstruct.secstruct_metric import SecstructMetric
from config import DATABASE_PATH
from util.gene_file_io import *

def get_features_filename(intron_class):
	return DATABASE_PATH + 'introns/' + intron_class + '/features.csv'

def get_zipper_stem_filename(intron_class):
	return DATABASE_PATH + 'introns/' + intron_class + '/zipper_stems.txt'

def build_intron_set(intron_class, feature_options={}, intron_options={}):
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

	add_bpp = False
	bpp_dir = ""
	if 'use_bpp' in feature_options.keys() and feature_options['use_bpp']:
		add_bpp = True
		bpp_dir = DATABASE_PATH + 'introns/' + intron_class + '/bpp_dir/'

	name_is_refseq = True
	get_ensembl_names = True
	if 'name_is_refseq' in intron_options.keys():
		name_is_refseq = intron_options['name_is_refseq']
	if 'get_ensembl_names' in intron_options.keys():
		get_ensembl_names = intron_options['get_ensembl_names']
	
	all_introns = IntronSet()
	all_introns.init_from_files(intron_seq_file, mfe_filename=mfe_file, \
		ens_filename=ens_file, name_is_refseq=name_is_refseq, \
		get_ensembl_names=get_ensembl_names, add_bpp=add_bpp, bpp_dir=bpp_dir)

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
		full_names_df = []
		for name in full_names:
			if name in features_df.columns:
				full_names_df += [name]
		features_df = features_df[full_names_df]

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

# Updates secondary structure files if they are not currently filled
# Raises an error if secondary structure flags are not set in the feature_options dictionary
def check_update_secstruct(intron_class, feature_options):
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


# Needs to also return all_introns because it can get updated to include
# the secondary structure here.
def get_feature_vals_update_secstruct(feature_name, intron_class, feature_options, all_introns):
	feature = get_feature_from_name(feature_name)
	
	# Make sure that intron set has secondary structures computed if needed
	if isinstance(feature, SecstructMetric):
		check_update_secstruct(intron_class, feature_options)
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
		if feature_vals != []:
			features_df[full_feature_name] = feature_vals
	
	features_file = get_features_filename(intron_class)	
	write_file = open(features_file, 'w')
	features_df.to_csv(path_or_buf=write_file, index=False)
	write_file.close()

# Get zipper stems for this intron class using the MFE
# If zipper stems have not been stored in the database cache, compute them fresh now
# and store to the database. Otherwise, retrieve them from the cache.
def get_zipper_stems(intron_class, feature_options):
	zipper_stem_filename = get_zipper_stem_filename(intron_class)

	force_eval = False
	if 'force_eval' in feature_options.keys():
		force_eval = feature_options['force_eval']
	verbose = False
	if 'verbose' in feature_options.keys():
		verbose = feature_options['verbose']

	if os.path.isfile(zipper_stem_filename) and not force_eval:
		zipper_stem_file_data = read_zipper_stem_file(zipper_stem_filename)
		return zipper_stem_file_data

	if ('secstruct_type' not in feature_options.keys()) or \
		(feature_options['secstruct_type'] != 'mfe'):
		raise RuntimeError("Feature options must include the secondary structure estimation type mfe")
	check_update_secstruct(intron_class, feature_options)
	all_introns = build_intron_set(intron_class, feature_options=feature_options, 
		intron_options={'get_ensembl_names': False})

	all_zipper_stems = []		
	for ii, intron in enumerate(all_introns.introns):
		if verbose:
			print("Feature: Zipper stems, Evaluating intron: %d of %d" % \
				(ii, len(all_introns.introns)))

		feature = get_feature_from_name("ZipperStemStartMetric")
		[has_dG, best_stem, _] = feature.has_stem_dG(intron.bp, intron.seq, intron.mfe)
		zipper_stem_data = process_zipper_stem_entry(has_dG, best_stem)
		all_zipper_stems += [zipper_stem_data]

	write_zipper_stem_file(zipper_stem_filename, all_zipper_stems)

	return all_zipper_stems

def get_feature_vals(feature, all_introns, feature_options):
	verbose = False
	if 'verbose' in feature_options.keys():
		verbose = feature_options['verbose']

	feature_vals = []
	for ii, intron in enumerate(all_introns.introns):
		if verbose:
			print("Feature: %s, Evaluating intron: %d of %d" % \
				(feature.name, ii, len(all_introns.introns)))
		try:
			feature_vals += [feature.apply(intron, feature_options)]
		except NotImplementedError:
			print("Feature: %s is not implemented" % feature.name)

	return feature_vals

# FeatureData holds options for features: secondary structure package, ens/mfe, print progress or no. 
# It is a dictionary from feature name to feature options
def get_features(feature_names, intron_class, feature_options_all={}, intron_options={}):
	all_introns = build_intron_set(intron_class, intron_options=intron_options)

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

