# Utilities for working with stockholm alignment files

import os

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]


# Check if sequence is empty (full of blanks), optionally in a 
# given range
def aln_seq_is_empty(sequence, aln_range=[]):
	if len(aln_range) == 0:
		aln_range = (0, len(sequence))

	is_empty = True
	for ii in range(aln_range[0], aln_range[1]):
		cur_char = sequence[ii]
		if not ((cur_char == ".") or (cur_char == "-")):
			is_empty = False
			break
	
	return is_empty


# Returns a list of intron positions in the form (5'SS, bp, 3'SS).
def get_introns(annotation_str):
	five_poses = find(annotation_str, '5')
	bp_poses = find(annotation_str, 'B')
	three_poses = find(annotation_str, '3')

	intron_list = []
	if (len(five_poses) == len(bp_poses)) and \
		(len(five_poses) == len(three_poses)):
		for ii in range(len(five_poses)):
			intron_list += \
				[(five_poses[ii]-2, bp_poses[ii]+3, three_poses[ii]+1)]
	return intron_list


# Remove gaps and update branchpoint position
def condense_seq(seq, pos_info):
	[five_ss, bp, three_ss] = pos_info
	full_seq = seq[five_ss:(three_ss+1)].upper()
	bp_pos = bp - five_ss
	new_bp_pos = -1
	new_seq = ""
	for ii in range(len(full_seq)):
		new_char = full_seq[ii]
		add_char = False
		if (new_char != '-') and (new_char != '.'):
			add_char = True
		if add_char:
			if ii == bp_pos:
				new_bp_pos = len(new_seq)
			new_seq += new_char
	return [new_seq, new_bp_pos]

# Remove gaps and update branchpoint position
def get_seq_pos_info(seq, pos_info):
	[five_ss, bp, three_ss] = pos_info
	full_seq = seq[five_ss:(three_ss+1)].upper()
	return full_seq

# For reading .stk file sequence lines across different \n file varieties. 
# Returns (number of lines skipped, starting position)
def get_skips(lines):
	num_skips = 0
	count = 0
	starting_pos = -1
	for line in lines:
		# This first one allows for up to 5 spaces in an "empty" line
		if len(line) < 5:
			num_skips += 1
		elif line[0] == '#':
			num_skips += 1
		elif line[0] == '/':
			num_skips += 1
		else:
			if (starting_pos == -1):
				starting_pos = count 
		count += 1
	return (num_skips, starting_pos)

def get_rnd_char():
	return random.sample(["A", "C", "G", "U"], 1)[0]

# Make a dictionary with:
# keys: species name 
# values: list of info for each intron: [(Gene name, intron seq, branchpoint pos, intron num), ...]
def make_species_seqs(alignment_dir):
	species_seqs = {}

	for filename in os.listdir(alignment_dir):
		if filename.endswith(".stk"):
			f = open(alignment_dir + filename)
			lines = f.readlines()
			f.close()

			gene_name = filename.split(".")[0]

			annotation_str = lines[len(lines)-2].split()[2]
			intron_pos_list = get_introns(annotation_str)
			print(filename)
			(num_skipped, starting_pos) = get_skips(lines)
			for ii in range(len(lines)-num_skipped):
				line = lines[ii + starting_pos]

				for jj in range(len(intron_pos_list)):
					intron_pos = intron_pos_list[jj]
					species_name = (line.split()[0]).split('_')[0]
					[species_seq, bp] = condense_seq(line.split()[1], intron_pos)

					if len(species_seq) == 0:
						continue

					new_entry = (gene_name, species_seq, bp, jj + 1)
					if species_name in species_seqs:
						species_seqs[species_name] += [new_entry]
					else:
						species_seqs[species_name] = [new_entry]

	return species_seqs


# Get stk alignment for this S. cer intron sequence
# Goes through all .stk alignment files in species_dir to find
# one with the sequence matching the intron sequence
# Returns a stockholm alignment dictionary with key: species name, value: sequence
# Includes scer species sequence with . and -'s
# If the intron sequence is not found in the .stk files, returns an empty dictionary
def get_intron_aln(intron_seq, ensembl_name, alignment_dir):
	intron_seq_rna = intron_seq.replace('T', 'U')

	aln_dict = {}

	# Assemble list of stockholm file alignments in the alignment directory, 
	# placing the most probable location of the intron's alignment first.
	all_filenames = []
	for filename in os.listdir(alignment_dir):
		if ensembl_name in filename:
			all_filenames = [filename] + all_filenames
		else:
			all_filenames += [filename]

	for filename in all_filenames:
		f = open(alignment_dir + filename)
		lines = f.readlines()
		f.close()

		annotation_str = lines[len(lines)-2].split()[2]
		intron_pos_list = get_introns(annotation_str)

		(num_skipped, starting_pos) = get_skips(lines)
		
		intron_hit = -1
		scer_hit = ""
		# First find out if the scer introns match the 
		# given intron_seq
		for ii in range(len(lines)-num_skipped):
			line = lines[ii + starting_pos]

			for jj in range(len(intron_pos_list)):
				intron_pos = intron_pos_list[jj]
				species_name = (line.split()[0]).split('_')[0]

				[species_seq, bp] = condense_seq(line.split()[1], intron_pos)
				
				if species_name == 'scer' and \
					intron_seq_rna == species_seq:
					intron_hit = jj
					scer_hit = get_seq_pos_info(line.split()[1], intron_pos)
		
		# If an scer intron matched the input intron_seq, 
		# assemble the dictionary of aligned sequences
		if intron_hit != -1:
			for ii in range(len(lines)-num_skipped):
				line = lines[ii + starting_pos]

				intron_pos = intron_pos_list[intron_hit]
				species_name = (line.split()[0]).split('_')[0]
				species_seq = get_seq_pos_info(line.split()[1], intron_pos)
				
				aln_dict[species_name] = species_seq.replace('U', 'T')

		# NOTE: this only keeps one paralog 
		if len(scer_hit) > 0:
			aln_dict['scer'] = scer_hit.replace('U', 'T')

		if intron_hit != -1:
			break

	return aln_dict

# Get index mapping from condensed intron sequence to 
# alignment sequence (alignment sequence includes . and - characters)
def get_idx_map_intron_to_aln_seq(intron_seq, aln_seq):
	idx_mapping = {}

	idx = 0
	for ii in range(len(aln_seq)):
		if (aln_seq[ii] == '-') or (aln_seq[ii] == '.'):
			continue
		if aln_seq[ii] != intron_seq[idx]:
			raise RuntimeError("Intron sequence does not match extended sequence from alignment")
		idx_mapping[idx] = ii
		idx += 1

	return idx_mapping
