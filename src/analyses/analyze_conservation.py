"""
For an intron set, plot histogram of number of orthologs and percent conservation across orthologs 

Use alignments from an alignment directory

Optionally also compute stats for zipper stems
"""
import argparse
import numpy as np
from matplotlib import pyplot as plt

from features_db import *
from aln_util import aln_seq_is_empty

def convert_intervals_to_aln(intron, intervals):
	new_intervals = []
	for interval in intervals:
		idx_mapping = intron.aln_idx_mapping
		new_start = idx_mapping[interval[0]]
		new_end = idx_mapping[interval[1]]
		new_intervals += [(new_start, new_end)]
	return new_intervals

def get_num_orthologs(intron, intervals):
	intervals = convert_intervals_to_aln(intron, intervals)

	num_orthologs = 0
	for species, species_seq in enumerate(intron.aln_dict.values()):
		is_empty = True
		for interval in intervals:
			if not aln_seq_is_empty(intron, interval):
				is_empty = False
				break
		if not is_empty:
			num_orthologs += 1
	return num_orthologs

# Get average conservation in range for all species with an ortholog for this intron
# NOTE: if a species has the full intron deleted, will not include in the calculation
def get_avg_conservation(intron, intervals, min_seqs=3)
	intervals = convert_intervals_to_aln(intron, intervals)

	total_conservation = 0

	all_seqs = []
	for species, species_seq in enumerate(intron.aln_dict.values()):
		if not aln_seq_is_empty(species_seq):
			all_seqs += [species_seq]

	if len(all_seqs) < min_seqs:
		return -1

	all_cons_vals = []
	for (start_idx, end_idx) in intervals:
		seq_arr = []
		for seq in all_seqs:
			arr += [list(seq[start_idx:end_idx])]
		seq_arr = np.array(arr)

		mask = (list(intron.seq) == seq_arr)
		cons_vals = np.sum(mask, axis=0)/len(all_seqs)
		all_cons_vals += [cons_vals]

	all_cons_vals = np.array(all_cons_vals)
	return np.mean(all_cons_vals)

# For now just provides the average conservation, but later could include things like
# number of intervals of at least __ length meeting a __% conservation cutoff
def get_stats(intron_class):
	all_introns = build_intron_set(intron_class)
	all_introns.fill_aln(alignment_dir)

	ortholog_stats = []
	cons_stats = []
	for intron in all_introns.introns:
		num_orthologs = get_num_orthologs(intron, [(0, len(intron.seq))])
		ortholog_stats += [num_orthologs]
		avg_conservation = get_avg_conservation(intron, [(0, len(intron.seq))])
		cons_stats += [avg_conservation]

	return [ortholog_stats, cons_stats]

def get_stats_zipper(intron_class):
	all_introns = build_intron_set(intron_class)
	all_introns.fill_aln(alignment_dir)

	feature_options = {'secstruct_pkg': 'Vienna', 
					'secstruct_type': 'mfe', 
					'verbose': True,
					'force_eval': False}

	zipper_stem_data = get_zipper_stems(intron_class, feature_options)

	ortholog_stats = []
	cons_stats = []

	for ii, intron in enumerate(all_introns.introns):
		zipper_stem = zipper_stem_data[ii]
		[seq1, seq2, _, _] = zipper_stem

		range1_start = intron.seq.index(seq1)
		range1_end = range1_start + len(intron.seq)
		range2_start = intron.seq.index(seq2)
		range2_end = range2_start + len(intron.seq)
		intervals = [(range1_start, range1_end), (range2_start, range2_end)]

		num_orthologs = get_num_orthologs(intron, intervals)
		avg_conservation = get_avg_conservation(intron, intervals)
		ortholog_stats += [num_orthologs]
		cons_stats += [avg_conservation]

	return [ortholog_stats, cons_stats]

def plot_ortholog_stats(ortholog_stats):
	plt.hist(ortholog_stats)
	plt.show()

def plot_conservation_stats(cons_stats):
	plt.hist(cons_stats)
	plt.show()

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Parameters for processing intron alignment data')
	parser.add_argument('alignment_dir', type=str, help='Path to directory storing alignments in Stockholm format')
	parser.add_argument('intron_class', type=str, help='Path to directory storing alignments in Stockholm format')
	parser.add_argument('--do_zipper', default=False, action='store_true', \
		 help='Do stats on zipper stems')	
	args = parser.parse_args()

	alignment_dir = args.alignment_dir
	intron_class = args.intron_class
	do_zipper = args.do_zipper


	[ortholog_stats, cons_stats] = get_stats(intron_class)
	plot_ortholog_stats(ortholog_stats)
	plot_conservation_stats(cons_stats)

	if do_zipper:
		[ortholog_stats, cons_stats] = get_stats_zipper(intron_class)
		plot_ortholog_stats(ortholog_stats)
		plot_conservation_stats(cons_stats)