from feature.secstruct.secstruct_metric import SecstructMetric
from core.secstruct import SecstructGraph

import networkx as nx

class MLDMetric(SecstructMetric):
	def __init__(self, max_ens = 10, name="MLDMetric"):
		self.max_ens = max_ens
		self.name = name

	def get_mld(self, secstruct):
		secstruct_graph = SecstructGraph(secstruct)
		G = secstruct_graph.G

		path_lens = nx.floyd_warshall(G, weight='weight')
		max_value = 0
		for cur_key in path_lens.keys():
			cur_max = max(path_lens[cur_key].values())
			if cur_max > max_value:
				max_value = cur_max
		return max_value/len(secstruct)

	def get_score_mfe(self, intron):
		start = intron.fivess_offset
		end = len(intron.seq) - intron.threess_offset
		mfe = intron.mfe[start:end]
		return self.get_mld(mfe)

	def get_score_ens(self, intron):
		start = intron.fivess_offset
		end = len(intron.seq) - intron.threess_offset

		total_mld = 0
		secstructs = intron.ens[:self.max_ens]
		for mfe in secstructs:
			trunc_mfe = mfe[start:end]
			total_mld += self.get_mld(trunc_mfe)

		return total_mld/len(secstructs)
