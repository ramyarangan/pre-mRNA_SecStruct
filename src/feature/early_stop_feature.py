from feature import Feature
from core.gene import STOP_CODONS
from core.gene import GeneSet

class EarlyStopFeature(Feature):
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

	def get_early_stops_and_failed_introns(intron, len_limit=-1):
		gene_set = GeneSet()
		genes_dict = gene_set.get_genes_dict()

		if intron.name not in genes_dict:
			raise RuntimeError("Intron not in ORF")

		gene = genes_dict[intron.name]
		if intron.seq not in gene.seq:
			raise RuntimeError("Intron not in gene sequence")

		original_gene_seq = gene.seq
		spliced_gene_seq = gene.seq.replace(intron.seq, '')
		
		stop1 = get_first_stop_in_frame(original_gene_seq)
		stop2 = get_first_stop_in_frame(spliced_gene_seq)
		
		if ((stop1 - len(intron.seq)) < stop2):
			late_stop = True

		if (stop2 < (stop1 - len(intron.seq))) and \
			((len_limit == -1) or (stop2 < len_limit)):
			early_stop = True

		threeprime_dist = len(spliced_gene_seq) - stop2

		return (early_stop, late_stop, threeprime_dist)

class HasEarlyStopFeature(Feature):
	def __init__(name="HasEarlyStopFeature"):
		self.name = name

	def apply(self, intron, feature_options):
		[is_early, is_late, threeprime_dist] = EarlyStop.get_early_stop_info(intron)
		if is_early:
			return 1
		return 0

class ThreeprimeDistStopFeature(Feature):
	def __init__(name="ThreeprimeDistStopFeature"):
		self.name = name

	def apply(self, intron, feature_options):
		[is_early, is_late, threeprime_dist] = EarlyStop.get_early_stop_info(intron)
		return threeprime_dist



