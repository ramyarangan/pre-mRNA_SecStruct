import networkx as nx

class Node():
	def __init__(self, base1, base2=-1, base_pair=False):
		self.base1 = base1
		self.base2 = base2
		self.base_pair = base_pair

class BigNode():
	def get_type(self):
		return 'base'

class Stem(BigNode):
	def __init__(self, strand1_nts, strand2_nts):
		if len(strand1_nts) != len(strand2_nts):
			raise RuntimeError("Stem needs to have equal # nts in strands")

		self.strand1_nts = strand1_nts
		self.strand2_nts = strand2_nts
		self.bpp = None

	def get_bpp(self, bpp_matrix):
		max_bpp = 0
		for ii, nt1 in enumerate(self.strand1_nts):
			nt2 = self.strand2_nts[ii]
			max_bpp = max(max_bpp, bpp_matrix[nt1][nt2])
		return max_bpp

	def set_bpp(self, bpp_matrix):
		self.bpp = self.get_bpp(bpp_matrix)

	def get_type(self):
		return 'stem'

	def len(self):
		return len(self.strand1_nts)

	def is_in(self, start_pos, end_pos):
		is_in = True
		for ii in self.strand1_nts:
			if ii < start_pos or ii > end_pos:
				is_in = False
		for ii in self.strand2_nts:
			if ii < start_pos or ii > end_pos:
				is_in = False
		return is_in

	def __str__(self):
		print_str = "Stem of length %d containing base-pairs:\n" % len(self.strand2_nts)
		for ii, n1 in enumerate(self.strand1_nts):
			print_str += "%d %d\n" % (n1, self.strand2_nts[ii])
		if self.bpp is not None:
			print_str += "Base pair probability: %f\n" % self.bpp
		return print_str

class Junction(BigNode):
	def __init__(self, nts):
		self.nts = nts

	def get_type(self):
		return 'junction'

class External(BigNode):
	def __init__(self, nts):
		self.nts = nts

	def get_type(self):
		return 'external'

"""
Secstruct class that stores the secondary structure as graph-based representations: 
type 1. each base-pair or single-stranded nt is a node
type 2. each junction, stem, or external ssRNA region is a node
"""
class SecstructGraph():
	def __init__(self, dotbracket):
		self.dotbracket = dotbracket
		stems, nt_to_stem_dict = self.get_stems()
		self.stems = stems
		self.nt_to_stem_dict = nt_to_stem_dict
		G, node_dict = self.get_graph()
		self.G = G
		self.node_dict = G
		self.stemG = self.get_stem_graph()

	def get_stems(self, stem_verbose=False):
		dotbracket = self.dotbracket
		stems = []

		# First build list of stems
		base1_list = []
		cur_stem = []
		nt_to_stem_dict = {}
		for ii, curchar in enumerate(dotbracket):
			in_stem = False
			if len(cur_stem) > 0 and curchar == ')':
				# We could be continuing a stem
				if base1_list[-1] == cur_stem[-1][0] - 1:
					in_stem = True

			if not in_stem and len(cur_stem) > 0:
				strand1_nts = [x[0] for x in cur_stem]
				strand2_nts = [x[1] for x in cur_stem]
				new_stem = Stem(strand1_nts, strand2_nts)
				if stem_verbose:
					print(new_stem)
				stems += [new_stem]
				for nt1 in strand1_nts:
					nt_to_stem_dict[nt1] = new_stem
				for nt2 in strand2_nts:
					nt_to_stem_dict[nt2] = new_stem
				cur_stem = []

			if curchar == '(':
				base1_list += [ii]

			if curchar == ')':
				cur_stem += [(base1_list[-1], ii)]
				base1_list = base1_list[:-1]

		if len(cur_stem) > 0:
			strand1_nts = [x[0] for x in cur_stem]
			strand2_nts = [x[1] for x in cur_stem]
			new_stem = Stem(strand1_nts, strand2_nts)
			stems += [new_stem]		
			for nt1 in strand1_nts:
				nt_to_stem_dict[nt1] = new_stem
			for nt2 in strand2_nts:
				nt_to_stem_dict[nt2] = new_stem

		return stems, nt_to_stem_dict

	# Node types: 
	#    Base-pair (base-pairing residues)
	#    Single nucleotide
	# Used for: MLD calculation (longest shortest paths)
	# Return type: networkx graph, node id to node object dict
	def get_graph(self, ss_weight=1, bp_weight=1):
		dotbracket = self.dotbracket
		base1_list = []
		node_dict = {}

		num_nodes = 0
		for ii, curchar in enumerate(dotbracket):
			if curchar == '.':
				num_nodes += 1
				node_dict[ii] = Node(ii)
			if curchar == '(':
				base1_list += [ii]
			if curchar == ')':
				base1 = base1_list[-1]
				base1_list = base1_list[:-1]
				new_node = Node(base1, base2=ii, base_pair=True)
				node_dict[ii] = new_node
				node_dict[base1] = new_node
				num_nodes += 1

		G = nx.Graph()

		for ii in range(len(dotbracket)):
			added_base = False
			if node_dict[ii].base_pair:
				# Base pair node
				if node_dict[ii].base1 == ii:
					G.add_node(ii)
					added_base = True
			else:
				# Loop node
				G.add_node(ii)
				added_base = True
			if added_base and ii > 0:
				prev_node = node_dict[ii - 1]
				weight = ss_weight
				if node_dict[ii].base_pair and prev_node.base_pair:
					weight = bp_weight
				G.add_edge(prev_node.base1, ii, weight=weight)

		return G, node_dict

	# Update graph with nodes for (cur_start, cur_end - 1) (inclusive)
	# nt_to_node_dict: for each nucleotide, the stem, junction, or external node it belongs to
	# dotbracket: dot bracket notation for the full structure
	# G: current graph
	# is_internal: False for a portion of the structure that is not in any stems 
	#              or internal to any stems
	def get_stem_graph_rec(self, cur_start, cur_end, nt_to_node_dict, \
		dotbracket, G, is_internal, verbose=False):
		# Add node for current portion of stem graph
		G.add_node(cur_start)

		if dotbracket[cur_start] == ')':
			raise RuntimeError("Unexpected recursion architecture: start index is a closing base-pair")

		# Case 1: stems. Here cur_start and cur_end are the ends of a stem; recursively update internal loops
		if dotbracket[cur_start] == '(':
			if verbose: 
				print("Doing stem at position: %d" % cur_start)
			if cur_start not in nt_to_node_dict.keys() or \
				nt_to_node_dict[cur_start].get_type() != 'stem':
				raise RuntimeError("Expected base-paired residue to be in a stem.")

			cur_stem = nt_to_node_dict[cur_start]
			
			stem_start = min(cur_stem.strand1_nts)
			junc_start = max(cur_stem.strand1_nts)
			junc_end = min(cur_stem.strand2_nts)
			stem_end = max(cur_stem.strand2_nts)
			
			if stem_start != cur_start:
				raise RuntimeError("Unexpected recursion architecture")

			# Update graph for nts between the 5' and 3' strands of stem
			# E.g. when Called on (((....)))... from (((((....)))...)), get node for ....
			node_id = self.get_stem_graph_rec(junc_start + 1, junc_end, nt_to_node_dict, \
				dotbracket, G, True, verbose=verbose)
			# Connect node for .... to ((()))
			G.add_edge(stem_start, node_id)
			
			node_id_2 = cur_start
			if cur_end > stem_end + 1:
				# This happens for internal loops with only 3' nucleotides or
				# external ssRNA outside of stems 
				# E.g. when Called on (((....)))... from (((((....)))...)), get node for ...
				node_id_2 = self.get_stem_graph_rec(stem_end + 1, cur_end, nt_to_node_dict, \
					dotbracket, G, is_internal, verbose=verbose)
				# Connect node for ... to ((()))
				G.add_edge(stem_start, node_id_2)

			if verbose:
				print("Finished stem at position: %d" % cur_start)

			if is_internal:
				# Called on (((....)))... from (((((....)))...))
				# Should return node for ... 
				# Called on (((....))) from ((...(((....)))...)) 
				# Should return node for ((()))
				return node_id_2
			# Called on (((....)))..... from ...(((....))).....
			# Should return node for ((()))
			return cur_start

		# all ssRNA in the loop / junctions / external nucleotides
		junc_nts = []

		# Case 2: external ssRNA
		if not is_internal:
			if verbose:
				print("Doing external ssRNA at position: %d" % cur_start)
			ssrna_end = cur_end
			cur_dotbracket = dotbracket[cur_start:cur_end]
			if '(' in cur_dotbracket: 
				ssrna_end = cur_dotbracket.index('(') + cur_start
				if ')' not in cur_dotbracket: 
					raise RuntimeError("Unexpected stem architecture")
				if dotbracket.index(')') < dotbracket.index('('):
					raise RuntimeError("Unexpected stem architectures")

			junc_nts = list(range(cur_start, ssrna_end))

			if cur_end > ssrna_end:
				# E.g. when called on ....(((...))) from ((....))....(((...)))
				node_id = self.get_stem_graph_rec(ssrna_end, cur_end, nt_to_node_dict, \
					dotbracket, G, False, verbose=verbose)
				# Connect .... to ((()))
				G.add_edge(cur_start, node_id)

			if verbose:
				print("Finished external ssRNA at position: %d" % cur_start)

			# Add new junction to datastructure
			new_external = External(junc_nts)
			for ii in junc_nts:
				nt_to_node_dict[ii] = new_external	

		# Case 3: internal ssRNA - loops and junctions
		else:
			if verbose:
				print("Doing internal loop / junction at position: %d" % cur_start)

			junc_nts = [] 
			neighbor_stems = [] # all stems internal to the junction
			cur_idx = cur_start
			while (cur_idx < cur_end):
				if dotbracket[cur_idx] == '.':
					junc_nts += [cur_idx]
					cur_idx += 1
				else:
					if cur_idx not in nt_to_node_dict.keys() or \
						nt_to_node_dict[cur_idx].get_type() != 'stem':
						raise RuntimeError("Expected base-paired residue to be in a stem.")
					neighor_stem = nt_to_node_dict[cur_idx]
					neighbor_stems += [neighor_stem]
					
					# Skip over all nucleotides in the new stem
					cur_idx = max(neighor_stem.strand2_nts) + 1

			# Recursion over all stems in the junction
			for neighbor_stem in neighbor_stems:
				start_idx = min(neighbor_stem.strand1_nts)
				end_idx = max(neighbor_stem.strand2_nts)
				# E.g. when called on ....((...))...(....).. from ((....((...))...(....)..))
				node_id = self.get_stem_graph_rec(start_idx, end_idx + 1, nt_to_node_dict, \
					dotbracket, G, True, verbose=verbose)
				# Connect ......... to (()) and ......... to ()
				G.add_edge(cur_start, node_id)

			if verbose:
				print("Finished internal loop / junction at position: %d" % cur_start)

			# Add new junction to datastructure
			new_junction = Junction(junc_nts)
			for ii in junc_nts:
				nt_to_node_dict[ii] = new_junction

		# When called on ....((...))...(....).. from ((....((...))...(....)..)) 
		# return start of .........
		# When called on ....(((...))) from ((....))....(((...)))
		# return start of ....
		return cur_start

	# Node types: 
	#      Stem (residues, length, bootstrapping probability) 
	#      Loop (all loop nts condensed into 1)
	# Used for: NWJ accounting, longest stem
	# Return type: networkx graph, dictionary connecting nucleotides to the stems/ junctions
	#              they belong in 
	def get_stem_graph(self, graph_verbose=False):
		dotbracket = self.dotbracket
		# Initiate stem graph building recursion
		G = nx.Graph()
		nt_to_node_dict = self.nt_to_stem_dict
		self.get_stem_graph_rec(0, len(dotbracket), nt_to_node_dict, dotbracket, G, \
			False, verbose=graph_verbose)
		return G

