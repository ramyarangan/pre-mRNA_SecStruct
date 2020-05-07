from feature.feature import Feature
from core.gene import STOP_CODONS
from core.gene import GeneSet
import numpy as np

class EarlyStop:
	# Returns the resulting ORF length based on location of first stop
	# Returns -1 if the STOP codon is not found
	def get_first_stop_in_frame(seq):
		ii = 0
		while ii < len(seq) - 2:
			next_codon = seq[ii:(ii+3)]
			if next_codon in STOP_CODONS:
				return (ii + 3)
			ii += 3
		return -1

	def get_early_stop_info(intron, genes_dict, len_limit=-1):
		if intron.name not in genes_dict:
			# Figure out why some gene annotations aren't right?
			# raise RuntimeError("Intron not in ORF: %s" % (intron.name))
			return (np.nan, np.nan, np.nan)

		gene = genes_dict[intron.name]
		if intron.seq not in gene.seq:
			raise RuntimeError("Intron not in gene sequence")

		original_gene_seq = gene.seq
		spliced_gene_seq = gene.seq.replace(intron.seq, '')
		
		stop1 = EarlyStop.get_first_stop_in_frame(original_gene_seq)
		stop2 = EarlyStop.get_first_stop_in_frame(spliced_gene_seq)
		
		early_stop = False
		late_stop = False

		if ((stop1 - len(intron.seq)) < stop2):
			late_stop = True

		if (stop2 < (stop1 - len(intron.seq))) and \
			((len_limit == -1) or (stop2 < len_limit)):
			early_stop = True

		threeprime_dist = len(spliced_gene_seq) - stop2
		
		return (early_stop, late_stop, threeprime_dist)

class HasEarlyStopFeature(Feature):
	def __init__(self, name="HasEarlyStopFeature"):
		self.name = name
		gene_set = GeneSet()
		genes_dict = gene_set.get_genes_dict()
		self.genes_dict = genes_dict

	def apply(self, intron, feature_options):
		[is_early, is_late, threeprime_dist] = \
			EarlyStop.get_early_stop_info(intron, self.genes_dict)
		if is_early:
			return 1
		return 0

class ThreeprimeDistStopFeature(Feature):
	def __init__(self, name="ThreeprimeDistStopFeature"):
		self.name = name
		gene_set = GeneSet()
		genes_dict = gene_set.get_genes_dict()
		self.genes_dict = genes_dict

	def apply(self, intron, feature_options):
		[is_early, is_late, threeprime_dist] = \
			EarlyStop.get_early_stop_info(intron, self.genes_dict)
		return threeprime_dist



