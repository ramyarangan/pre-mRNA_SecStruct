from feature.secstruct.secstruct_metric import SecstructMetric

class LocalizationMetric(SecstructMetric):
	# Scores 2D graph localization between pair1 and pair2
	def __init__(self, start_bp=False, end_bp=False, name="LocalizationMetric"):
		self.INTMAX = 100000
		self.start_bp = start_bp
		self.end_bp = end_bp
		self.name = name

	def build_nodes(self, secstruct):
		node_to_ids = {}
		ids_to_node = {}
		cur_node = 0
		bp_stack = []
		for ii in range(0, len(secstruct)):
			if (secstruct[ii] == '.'):
				node_to_ids[cur_node] = [ii]
				ids_to_node[ii] = cur_node
				cur_node = cur_node + 1
			if (secstruct[ii] == "("):
				bp_stack = bp_stack + [ii]
			if (secstruct[ii] == ")"):
				bp_start = bp_stack[-1]
				bp_stack = bp_stack[0:-1]
				node_to_ids[cur_node] = [bp_start, ii]
				ids_to_node[bp_start] = cur_node
				ids_to_node[ii] = cur_node
				cur_node = cur_node + 1
		return (node_to_ids, ids_to_node)

	def build_edges(self, ids_to_node, seq_length):
		edge_map = {} # node: list of adj nodes
		for ii in range(0, seq_length - 1):
			node1 = ids_to_node[ii]
			node2 = ids_to_node[ii + 1]
			if (node1 != node2):
				if node1 in edge_map:
					edge_map[node1] = list(set(edge_map[node1] + [node2]))
				else:
					edge_map[node1] = [node2]
				if node2 in edge_map:
					edge_map[node2] = list(set(edge_map[node2] + [node1]))
				else:
					edge_map[node2] = [node1]
		return edge_map

	def shortest_path(self, ids_to_node, edge_map, id_1, id_2, num_nodes):
		dists = [self.INTMAX] * num_nodes
		queue = [self.INTMAX] * num_nodes
		num_in_queue = 1
		cur_node = ids_to_node[id_1]
		dists[cur_node] = 0
		queue[cur_node] = 0
		target_idx = ids_to_node[id_2]
		final_dist = -1
		while (num_in_queue > 0):
			# Get min element from queue
			# Ignore elements that have been removed
			min_idx = -1
			min_val = self.INTMAX
			for ii in range(0, num_nodes):
				if (queue[ii] < min_val) and (queue[ii] != -1):
					min_idx = ii
					min_val = queue[ii]
			
			# Remove a value from the queue
			num_in_queue = num_in_queue - 1
			queue[min_idx] = -1
			if (min_idx == target_idx):
				final_dist = min_val
				break

			# Go through neighbors and update their values in queue
			neighbors = edge_map[min_idx]
			for neighbor in neighbors:
				new_dist = min_val + 1
				if (queue[neighbor] == self.INTMAX):
					num_in_queue = num_in_queue + 1
				if (new_dist < queue[neighbor]):
					queue[neighbor] = new_dist
		return final_dist

	def get_mean_dist(self, secstructs, start, end):
		bp_dists = []
		for secstruct in secstructs:
			(node_to_ids, ids_to_node) = self.build_nodes(secstruct)
			edge_map = self.build_edges(ids_to_node, len(secstruct))
			num_nodes = len(node_to_ids.keys())
			bp_dist = self.shortest_path(ids_to_node, edge_map, start, end, num_nodes)
			bp_dists = bp_dists + [bp_dist]
		return sum(bp_dists)/len(secstructs)

	def get_score_mfe(self, intron):
		start = intron.fivess_offset
		end = len(intron.seq) - intron.threess_offset
		if self.start_bp:
			start = intron.bp - intron.fivess_offset
		if self.end_bp:
			end = intron.bp - intron.fivess_offset
		return self.get_mean_dist([intron.mfe], start, end)

	def get_score_ens(self, intron):
		start = intron.fivess_offset
		end = len(intron.seq) - intron.threess_offset
		if self.start_bp:
			start = intron.bp - intron.fivess_offset
		if self.end_bp:
			end = intron.bp - intron.fivess_offset
		return self.get_mean_dist(intron.ens, start, end)
