from feature.secstruct.secstruct_metric import SecstructMetric
from core.secstruct import SecstructGraph

import networkx as nx

class MLDMetric(SecstructMetric):
	def __init__(self, max_ens = 10, name="MLDMetric"):
		self.max_ens = max_ens
		self.name = name

	def get_mld(self, secstruct, fivess_offset, threess_offset):
		secstruct_graph = SecstructGraph(secstruct)
		G = secstruct_graph.G
		node_dict = secstruct_graph.node_dict

		path_lens = nx.floyd_warshall(G, weight='weight')
		max_value = 0

		for cur_key in path_lens.keys():
			end_idx = len(secstruct) - threess_offset
			if cur_key < fivess_offset or cur_key > end_idx:
				continue
			start_node = node_dict[fivess_offset].base1
			end_node = node_dict[end_idx - 1].base1
			cur_max = min(path_lens[cur_key][start_node], \
				path_lens[cur_key][end_node])
			max_value = max(cur_max, max_value)
		#	cur_max = max(path_lens[cur_key].values())
		#	if cur_max > max_value:
		#		max_value = cur_max
		return max_value/len(secstruct)


	def get_score_mfe(self, intron):
		return self.get_mld(intron.mfe, intron.fivess_offset, intron.threess_offset)

	def get_score_ens(self, intron):
		total_mld = 0
		secstructs = intron.ens[:self.max_ens]
		for mfe in secstructs:
			total_mld += self.get_mld(mfe, intron.fivess_offset, intron.threess_offset)

		return total_mld/len(secstructs)
