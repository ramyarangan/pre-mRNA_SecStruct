"""
For an intron set, plot histogram of number of orthologs and percent conservation across orthologs 

Use alignments from an alignment directory

Optionally also compute stats for zipper stems

Example usage: 
python analyses/analyze_conservation.py ../database/alignments/hooks_alignments/ standard_allsize_min_50_max_600  --do_zipper
"""
import os
import argparse
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster

from util.features_db import *
from util.aln_util import aln_seq_is_empty

# In order of evolutionary distance
SPECIES_NAMES = ['scer', 'smik', 'skud', 'suva', 'cgla', 'kafr', 'knag', 'ncas', \
	'ndai', 'tbla', 'tpha', 'kpol', 'zrou', 'tdel', 'klac', 'agos', 'ecym', 'sklu', \
	'kthe', 'kwal']
# use cutoff of 40 zipper stems
# SPECIES_NAMES = ['scer', 'smik', 'skud', 'suva', 'kafr', 'knag', 'ncas', 'ndai', \
# 	'kpol', 'klac', 'ecym', 'sklu', 'kthe', 'kwal']

def convert_intervals_to_aln(intron, intervals):
	new_intervals = []
	for interval in intervals:
		idx_mapping = intron.aln_idx_mapping
		new_start = idx_mapping[interval[0]]
		new_end = idx_mapping[interval[1]-1]
		new_intervals += [(new_start, new_end+1)]
	return new_intervals

def get_num_orthologs(intron, intervals):
	intervals = convert_intervals_to_aln(intron, intervals)

	num_orthologs = 0
	for species, species_seq in enumerate(intron.aln_dict.values()):
		is_empty = True
		for interval in intervals:
			if not aln_seq_is_empty(species_seq, aln_range=interval):
				is_empty = False
				break
		if not is_empty:
			num_orthologs += 1
	return num_orthologs

# Get average conservation in range for all species with an ortholog for this intron
# NOTE: if a species has the full intron deleted, will not include in the calculation
def get_avg_conservation(intron, intervals, min_seqs=3):
	intervals = convert_intervals_to_aln(intron, intervals)

	total_conservation = 0

	all_seqs = []
	scer_seq = ""
	for species, species_seq in intron.aln_dict.items():
		if not aln_seq_is_empty(species_seq):
			all_seqs += [species_seq]
			if species == 'scer':
				scer_seq = species_seq

	if len(all_seqs) < min_seqs:
		return -1

	all_cons_vals = np.array([])
	for (start_idx, end_idx) in intervals:
		seq_arr = []
		for seq in all_seqs:
			seq_arr += [list(seq[start_idx:end_idx])]
		seq_arr = np.array(seq_arr)

		scer_int = scer_seq[start_idx:end_idx]
		mask = (list(scer_int) == seq_arr)
		cons_vals = np.sum(mask, axis=0)/len(all_seqs)
		all_cons_vals = np.append(all_cons_vals, cons_vals)

	return np.mean(all_cons_vals)

def has_alignments(intron):
	if len(intron.aln_idx_mapping.keys()) == 0:
		return False
	return True

# For now just provides the average conservation, but later could include things like
# number of intervals of at least __ length meeting a __% conservation cutoff
def get_stats(all_introns, intron_class):
	ortholog_stats = []
	cons_stats = []
	for intron in all_introns.introns:
		if not has_alignments(intron):
			print("Skipping intron without alignment: %s" % intron.print_string())
			continue
		num_orthologs = get_num_orthologs(intron, [(0, len(intron.seq))])
		ortholog_stats += [num_orthologs]
		avg_conservation = get_avg_conservation(intron, [(0, len(intron.seq))])
		cons_stats += [avg_conservation]

	return [ortholog_stats, cons_stats]

def get_stats_zipper(all_introns, intron_class):
	feature_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'mfe', 
					'verbose': True,
					'force_eval': False}

	zipper_stem_data = get_zipper_stems(intron_class, feature_options)

	ortholog_stats = []
	cons_stats = []

	for ii, intron in enumerate(all_introns.introns):
		if not has_alignments(intron):
			print("Skipping intron without alignment: %s" % intron.print_string())
			continue		
		zipper_stem = zipper_stem_data[ii]
		[seq1, seq2, _, _] = zipper_stem

		if len(seq1) > 0 and len(seq2) > 0:
			range1_start = intron.seq.index(seq1)
			range1_end = range1_start + len(seq1)
			range2_start = intron.seq.index(seq2)
			range2_end = range2_start + len(seq2)
			intervals = [(range1_start, range1_end), (range2_start, range2_end)]

			num_orthologs = get_num_orthologs(intron, intervals)
			avg_conservation = get_avg_conservation(intron, intervals)
			ortholog_stats += [num_orthologs]
			cons_stats += [avg_conservation]

	return [ortholog_stats, cons_stats]

def plot_ortholog_stats(ortholog_stats):
	counts, bins, bars = plt.hist(ortholog_stats, range=[0, 20])
	plt.show()
	print(counts)
	print(bins)

def plot_conservation_stats(cons_stats):
	counts, bins, bars = plt.hist(cons_stats, range=[0, 1])
	plt.show()
	print(counts)
	print(bins)

def plot_heatmap(data_arr, data_mask, labels, num_homolog_bars, perc_bars, colormap='YlGnBu_r'):
	row_linkage = hierarchy.linkage(distance.pdist(data_arr))
	permutation = hierarchy.leaves_list(row_linkage)
	data_arr[:] = data_arr[permutation,:]
	data_mask[:] = data_mask[permutation,:]
	num_homolog_bars = np.array(num_homolog_bars)
	num_homolog_bars[:] = num_homolog_bars[permutation]
	num_homolog_bars = list(num_homolog_bars)

	fig, ((cbar_ax, ax1, hbar_ax), (dummy_ax, ax2, dummy_ax2)) = \
		plt.subplots(nrows=2, ncols=3, 
		figsize=(4, 8), sharex='col',
		gridspec_kw={'height_ratios': [10, 1], 'width_ratios': [1, 10, 2]})

	colormap = plt.cm.get_cmap(colormap)
	colormap.set_bad('gray')
	g = sns.heatmap(data_arr, mask=data_mask, cmap=colormap, xticklabels=False, \
		yticklabels=False, ax=ax1, cbar_ax=cbar_ax)

	ax2.set_xlabel("Species")
	x_tick_pos = [i + 0.5 for i in range(len(labels))]
	ax2.bar(x_tick_pos, perc_bars, align='center', color='grey')
	ax2.set_ylim(0, 1)
	ax2.set_xticks(x_tick_pos)
	ax2.set_xticklabels(labels, Rotation=45, Fontsize=8)

	y_tick_pos = [i + 0.5 for i in range(len(num_homolog_bars))]
	hbar_ax.barh(y_tick_pos, num_homolog_bars[::-1], align='center', color='grey')
	hbar_ax.set_ylim(0, len(num_homolog_bars))
	hbar_ax.set_yticks([])

	dummy_ax.axis('off')
	dummy_ax2.axis('off')

	plt.tight_layout()
	plt.show()

def plot_paralogs_separately(dG_species_dict, all_gene_names, \
	do_binary=False, dG_cutoff=0):

	# Get heatmap data - each aligned intron can match up to two introns/species 
	# two introns whenever this gene was duplicated in the WGD event
	dG_arr = np.zeros((len(all_gene_names), len(SPECIES_NAMES) * 2))
	dG_mask = np.full((len(all_gene_names), len(SPECIES_NAMES) * 2), False)
	perc_bars = []
	num_homolog_bars = [0] * len(all_gene_names)
	for jj, species in enumerate(SPECIES_NAMES):
		species_zipper_cnt = 0
		species_total = 0
		for ii, name in enumerate(all_gene_names):
			if name not in dG_species_dict[species]:
				dG_mask[ii, jj] = True
				dG_mask[ii, jj + len(SPECIES_NAMES)] = True
				continue
			dG_vals = dG_species_dict[species][name]
			if len(dG_vals) > 2:
				raise RuntimeError("Do not expect more than one paralog per gene")
			dG_arr[ii, jj] = dG_vals[0]
			if do_binary:
				dG_arr[ii, jj] = (float(dG_arr[ii, jj]) < dG_cutoff)
			if len(dG_vals) > 1:
				dG_arr[ii, jj + len(SPECIES_NAMES)] = dG_vals[1]
				if do_binary:
					dG_arr[ii, jj + len(SPECIES_NAMES)] = \
						(float(dG_vals[1]) < dG_cutoff)
			else:
				dG_mask[ii, jj + len(SPECIES_NAMES)] = True
			if do_binary: 
				species_zipper_cnt += dG_arr[ii, jj] + \
					dG_arr[ii, jj + len(SPECIES_NAMES)]
				if dG_arr[ii, jj] or (dG_arr[ii, jj + len(SPECIES_NAMES)]):
					num_homolog_bars[ii] += 1
			else:
				species_zipper_cnt += (dG_arr[ii, jj] < 0)
				species_zipper_cnt += (dG_arr[ii, jj + len(SPECIES_NAMES)] < 0)
				if (dG_arr[ii, jj] < 0) or \
					(dG_arr[ii, jj + len(SPECIES_NAMES)] < 0):
					num_homolog_bars[ii] += 1
			if not dG_mask[ii, jj]:
				species_total += 1
			if not dG_mask[ii, jj + len(SPECIES_NAMES)]:
				species_total += 1
		perc_bars += [species_zipper_cnt/species_total]
	labels = [0] * len(SPECIES_NAMES) * 2
	for ii, species in enumerate(SPECIES_NAMES):
		labels[ii] = species
		labels[ii + len(SPECIES_NAMES)] = species

	colormap = 'YlGnBu_r'
	if do_binary:
		colormap='YlGnBu'
	plot_heatmap(dG_arr, dG_mask, labels, num_homolog_bars, perc_bars, colormap=colormap)

def plot_best_zipper_paralog(dG_species_dict, all_gene_names, \
	do_binary=False, dG_cutoff=0):

	# Get heatmap data - each aligned intron can match up to two introns/species 
	# two introns whenever this gene was duplicated in the WGD event
	dG_arr = np.zeros((len(all_gene_names), len(SPECIES_NAMES)))
	dG_mask = np.full((len(all_gene_names), len(SPECIES_NAMES)), False)
	perc_bars = []
	num_homolog_bars = [0] * len(all_gene_names)
	for jj, species in enumerate(SPECIES_NAMES):
		species_zipper_cnt = 0
		total_intron_cnt = 0
		for ii, name in enumerate(all_gene_names):
			if name not in dG_species_dict[species]:
				dG_mask[ii, jj] = True
				continue
			dG_vals = dG_species_dict[species][name]
			if len(dG_vals) > 2:
				raise RuntimeError("Do not expect more than one paralog per gene")
			dG_arr[ii, jj] = dG_vals[0]
			if len(dG_vals) > 1:
				dG_arr[ii, jj] = min(dG_vals[0], dG_vals[1])
			if do_binary:
				dG_arr[ii, jj] = (float(dG_arr[ii, jj]) < dG_cutoff)
			if not do_binary and dG_arr[ii, jj] < 0: 
				species_zipper_cnt += 1
				num_homolog_bars[ii] += 1
			if do_binary and dG_arr[ii, jj]: 
				species_zipper_cnt += 1
				num_homolog_bars[ii] += 1
			total_intron_cnt += 1				

		species_zipper_perc = species_zipper_cnt/total_intron_cnt
		print("Zipper stems for %s: %f, %d" % (species, \
			species_zipper_perc, species_zipper_cnt))
		perc_bars += [species_zipper_perc]

	colormap = 'YlGnBu_r'
	if do_binary:
		colormap='Greens'
	plot_heatmap(dG_arr, dG_mask, SPECIES_NAMES, num_homolog_bars, perc_bars, colormap=colormap)

def plot_histogram_num_zipperstems(dG_species_dict, all_gene_names, dG_cutoff=0):
	num_zipper_stems = []
	for gene_name in all_gene_names:
		cnt_zipper_stem = 0
		for species in SPECIES_NAMES:
			if gene_name not in dG_species_dict[species].keys():
				continue
			dG_vals = dG_species_dict[species][gene_name]
			if len(dG_vals) > 0:
				if float(dG_vals[0]) < dG_cutoff:
					cnt_zipper_stem += 1
			if len(dG_vals) > 1:
				if float(dG_vals[1]) < dG_cutoff:
					cnt_zipper_stem += 1
		num_zipper_stems += [cnt_zipper_stem]

	bins = np.arange(0, max(num_zipper_stems))
	ticks = np.arange(0, max(num_zipper_stems), 2)
	counts, bins, bars = plt.hist(num_zipper_stems, \
		color='blue', bins=bins, alpha=0.5, rwidth=0.85)

	plt.xticks(ticks=ticks, labels=ticks)
	plt.xlabel("Number of homologous introns with zipper stems")
	plt.ylabel("Frequency")
	plt.title("Distribution of zipper stems in introns in Saccharomyces genus")
	plt.show()

	print(counts)
	print(bins)

def heatmap_zipper_stem_data(intron_class):
	intron_class = 'species_hooks_genes_named/' + intron_class
	feature_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'mfe', 
					'verbose': True,
					'force_eval': False}

	intron_path = os.path.join(DATABASE_PATH, 'introns/' + intron_class)

	species_names = []
	for species_name in os.listdir(intron_path):
		species_names += [species_name]

	# Retrieve zipper stem data and names
	dG_species_dict = {}
	gene_name_set = set() # Superset of all introns across all species - note that 
					# some species may not have a representative for an intron
	for species_name in species_names:
		dG_species_dict[species_name] = {}
		intron_class_species = intron_class + '/' + species_name
		zipper_stem_data = get_zipper_stems(intron_class_species, feature_options)
		print(zipper_stem_data)
		dGs = [x[3] for x in zipper_stem_data]
		all_introns = build_intron_set(intron_class_species, \
			intron_options={'name_is_refseq': False})
		names = [intron.name for intron in all_introns.introns]
		for ii, name in enumerate(names):
			if name in dG_species_dict[species_name].keys():
				dG_species_dict[species_name][name] += [dGs[ii]]
			else:
				dG_species_dict[species_name][name] = [dGs[ii]]
		gene_name_set = gene_name_set.union(set(names))
	all_gene_names = list(gene_name_set)

	plot_histogram_num_zipperstems(dG_species_dict, all_gene_names)
	# plot_paralogs_separately(dG_species_dict, all_gene_names)
	# plot_paralogs_separately(dG_species_dict, all_gene_names, do_binary=True)
	plot_best_zipper_paralog(dG_species_dict, all_gene_names)
	plot_best_zipper_paralog(dG_species_dict, all_gene_names, do_binary=True)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Parameters for processing intron alignment data')
	parser.add_argument('alignment_dir', type=str, help='Path to directory storing alignments in Stockholm format')
	parser.add_argument('intron_class', type=str, help='Intron class to analyze conservation values for')
	parser.add_argument('--do_zipper', default=False, action='store_true', \
		 help='Do stats on zipper stems')
	parser.add_argument('--do_stats', default=False, action='store_true', \
		 help='Do stats on alignment')		
	parser.add_argument('--do_zipper_species', default=False, action='store_true', \
		 help='Get zipper stems for all species')	
	args = parser.parse_args()

	alignment_dir = args.alignment_dir
	intron_class = args.intron_class
	do_zipper = args.do_zipper
	do_stats = args.do_stats
	do_zipper_species= args.do_zipper_species

	if do_stats or do_zipper:
		all_introns = build_intron_set(intron_class)
		all_introns.fill_aln(alignment_dir)

		if do_stats:
			[ortholog_stats, cons_stats] = get_stats(all_introns, intron_class)
			plot_ortholog_stats(ortholog_stats)
			plot_conservation_stats(cons_stats)

		if do_zipper:
			[ortholog_stats, cons_stats] = get_stats_zipper(all_introns, intron_class)
			plot_ortholog_stats(ortholog_stats)
			plot_conservation_stats(cons_stats)

	if do_zipper_species:
		heatmap_zipper_stem_data(intron_class)


