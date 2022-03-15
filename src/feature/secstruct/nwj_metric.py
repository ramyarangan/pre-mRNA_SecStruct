from feature.secstruct.secstruct_metric import SecstructMetric
from core.secstruct import SecstructGraph

import networkx as nx

class NWJMetric(SecstructMetric):
	def __init__(self, max_ens = 1000, name="NWJMetric"):
		self.max_ens = max_ens
		self.name = name

	def get_num_nwj(self, secstruct, loop_cutoff=20, min_n=3):
		secstruct_graph = SecstructGraph(secstruct)
		stemG = secstruct_graph.stemG
		nt_to_node_dict = secstruct_graph.nt_to_stem_dict

		nwjs = []

		for nt in nx.nodes(stemG):
			if nt not in nt_to_node_dict.keys():
				raise RuntimeError("Unexpected missing BigNode for nt: %s" % nt)

			if nt_to_node_dict[nt].get_type() != 'junction':
				continue

			# Check that there are at least min_n neighbors
			if len(list(nx.neighbors(stemG, nt))) != min_n:
				continue

			# Check that loop is at most loop_cutoff length
			if len(nt_to_node_dict[nt].nts) > loop_cutoff:
				continue

			nwjs += [nt_to_node_dict[nt]]

		return len(nwjs)

	def get_score_mfe(self, intron):
		return self.get_num_nwj(intron.mfe)

	def get_score_ens(self, intron):
		total_nwj = 0
		secstructs = intron.ens[:self.max_ens]
		for mfe in secstructs:
			total_nwj += self.get_num_nwj(mfe)

		return total_nwj/len(secstructs)
