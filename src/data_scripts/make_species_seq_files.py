import sys
import os
import random 
import numpy as np

alignment_dir = sys.argv[1]
seq_dir = sys.argv[2]
control_dir = ""
if len(sys.argv) > 3:
	control_dir = sys.argv[3]

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]


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


def get_seq(seq, pos_info):
	[five_ss, bp, three_ss] = pos_info
	full_seq = seq[five_ss:(three_ss+1)].upper()
	return full_seq


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


def make_mutations(scer_intron, scer_intron_exp, species_seq):
	# Tabulate modification counts
	insertions = 0
	deletions = 0
	mutations = 0
	for ii in range(len(scer_intron_exp)):
		if (scer_intron_exp[ii] == ".") or (scer_intron_exp[ii] == "-"):
			if (species_seq[ii] != ".") and (species_seq[ii] != "-"):
				insertions += 1
		elif (species_seq[ii] == ".") or (species_seq[ii] == "-"):
			if (scer_intron_exp[ii] != ".") or (scer_intron_exp[ii] != "-"):
				deletions += 1
		elif (scer_intron_exp[ii] != species_seq[ii]):
			mutations += 1
	print(str(insertions) + " " + str(deletions) + " " + str(mutations))
	# Make a np array of characters
	intron = np.array(list(scer_intron))
	
	# Do the deletions
	keep_locs = random.sample(range(len(scer_intron)), len(scer_intron) - deletions)
	keep_locs.sort()
	intron = intron[keep_locs]
	
	# Make space for the insertions
	insert_intron = np.array([" "] * (len(intron) + insertions))
	intron_locs = random.sample(range(len(insert_intron)), len(intron))
	intron_locs.sort()
	insert_intron[intron_locs] = intron
	intron = insert_intron

	# Do the mutations and insertions
	for ii in range(len(intron)):
		if intron[ii] == " ":
			intron[ii] = get_rnd_char()

	mut_locs = random.sample(range(len(intron)), mutations)
	for mut_loc in mut_locs:
		old_char = intron[mut_loc]
		while old_char == intron[mut_loc]:
			intron[mut_loc] = get_rnd_char()

	# Return in string form
	return "".join(intron)



# First make species: [(Gene name, intron seq, branchpoint pos, intron num), ...] dictionary.
def make_control_seqs(alignment_dir):
	control_seqs = {}

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

			for ii in range(len(intron_pos_list)):
				intron_pos = intron_pos_list[ii]

				scer_intron_exp = ""
				scer_intron = "" 
				for jj in range(len(lines)-num_skipped):
					line = lines[jj + starting_pos]
					species_name = (line.split()[0]).split('_')[0]
					if (species_name == "scer"):
						scer_intron_exp = get_seq(line.split()[1], intron_pos)
						[scer_intron, _] = condense_seq(line.split()[1], intron_pos)

				for jj in range(len(lines)-num_skipped):
					line = lines[jj + starting_pos]
					species_name = (line.split()[0]).split('_')[0]
					species_seq = get_seq(line.split()[1], intron_pos)
					[condensed_seq, bp] = condense_seq(line.split()[1], intron_pos)
					print(species_name)
					if len(condensed_seq) == 0:
						continue
					mutated_seq = make_mutations(scer_intron, scer_intron_exp, species_seq)
					# print(len(mutated_seq)) # Useful check that mutations are reasonable
					# print(len(condensed_seq))

					new_entry = (gene_name, mutated_seq, bp, ii + 1)
					if species_name in control_seqs:
						control_seqs[species_name] += [new_entry]
					else:
						control_seqs[species_name] = [new_entry]

	return control_seqs


# First make species: [(Gene name, intron seq, branchpoint pos, intron num), ...] dictionary.
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


def write_seqdat_files(species_seqs, seq_dir):
	# Get number of introns for S. cer for labeling. 
	intron_cnts = {}
	for intron in species_seqs['scer']:
		gene_name = intron[0]
		if gene_name in intron_cnts:
			intron_cnts[gene_name] = max(intron_cnts[gene_name], intron[3])
		else:
			intron_cnts[gene_name] = 1

	for species_name, introns in species_seqs.items():
		f = open(seq_dir + species_name + '.dat', 'w')

		for intron in introns:
			gene_name = intron[0]
			if gene_name not in intron_cnts:
				continue

			if intron_cnts[gene_name] > 1:
				gene_name += '_' + str(intron[3])

			info_str = str(intron[2]) + ' -1 -1 -1 ' + gene_name
			f.write('%s\n' % info_str)
			f.write('%s\n' % intron[1])

		f.close()


species_seqs = make_species_seqs(alignment_dir)
write_seqdat_files(species_seqs, seq_dir)
if (control_dir != ""):
	control_seqs = make_control_seqs(alignment_dir)
	write_seqdat_files(control_seqs, control_dir)