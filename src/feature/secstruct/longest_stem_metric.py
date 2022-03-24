from feature.secstruct.secstruct_metric import SecstuctMetricBPP
from core.secstruct import SecstructGraph

import networkx as nx

class LongestStemMetric(SecstuctMetricBPP):
	def __init__(self, max_ens = 1000, name="LongestStemMetric"):
		self.max_ens = max_ens
		self.name = name
	
	def get_longest_stem(self, secstruct, loop_cutoff=10, do_bpp=False, \
		bpp_thresh=0.7, bpp_len_cutoff=4, bpp_matrix=None):
		secstruct_graph = SecstructGraph(secstruct)
		stemG = secstruct_graph.stemG
		nt_to_node_dict = secstruct_graph.nt_to_stem_dict

		# All junctions must be 2WJ
		# < 10 loop nucleotides
		tested_stems = set()
		longest_stem = []
		longest_len = 0
		for nt in nx.nodes(stemG):
			# Only explore each node once
			if nt in tested_stems: 
				continue
			tested_stems.add(nt)

			if nt not in nt_to_node_dict.keys():
				raise RuntimeError("Unexpected missing BigNode for nt: %s" % nt)

			# Extend current stem as far as possible
			cur_longest_stem = []
			cur_longest_len = 0

			if nt_to_node_dict[nt].get_type() != 'stem':
				continue
			
			cur_stem = nt_to_node_dict[nt]

			cur_longest_stem += [cur_stem]
			cur_longest_len += cur_stem.len()

			passes_bpp = True
			if do_bpp:
				passes_bpp = False
				if cur_stem.get_bpp(bpp_matrix) < bpp_thresh:
					continue
				if cur_stem.len() > bpp_len_cutoff:
					passes_bpp = True

			# Extend stem as far as possible
			neighbors = list(nx.neighbors(stemG, nt))
			while len(neighbors) > 0:
				cur_neighbor = neighbors[0]
				neighbors = neighbors[1:]

				if cur_neighbor in tested_stems:
					continue
				
				tested_stems.add(cur_neighbor)

				neighbor_node = nt_to_node_dict[cur_neighbor]

				if neighbor_node.get_type() == 'stem':
					if do_bpp:
						if neighbor_node.get_bpp(bpp_matrix) < bpp_thresh:
							continue
						if neighbor_node.len() > bpp_len_cutoff:
							passes_bpp = True
					cur_longest_stem += [neighbor_node]
					cur_longest_len += neighbor_node.len()

				if neighbor_node.get_type() == 'loop':
					if len(neighbor_node.nts) > loop_cutoff:
						continue
					next_neighbors = nx.neighbors(stemG, cur_neighbor)
					if len(next_neighbors) > 2:
						continue

				if neighbor_node.get_type() == 'external':
					continue
				# Add all neighbors to the neighbor list to search for
				# extensions to the stem. 
				new_neighbors = nx.neighbors(stemG, cur_neighbor)
				for candidate in new_neighbors:
					if candidate not in tested_stems:
						neighbors += [candidate]

			if passes_bpp and cur_longest_len > longest_len:
				longest_len = cur_longest_len
				longest_stem = cur_longest_stem

		return longest_stem, longest_len

	def get_score_mfe(self, intron):
		_, longest_len = self.get_longest_stem(intron.mfe)
		return longest_len

	def get_score_mfe_bpp(self, intron):
		_, longest_len = self.get_longest_stem(intron.mfe, do_bpp=True, \
			bpp_matrix=intron.bpp)
		return longest_len

	def get_score_ens(self, intron):
		total_longest_len = 0
		secstructs = intron.ens[:self.max_ens]
		for mfe in secstructs:
			_, longest_len = self.get_longest_stem(mfe)
			total_longest_len += longest_len

		return total_longest_len/len(secstructs)
