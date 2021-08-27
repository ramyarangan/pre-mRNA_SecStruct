import sys
import math
import random
import numpy as np 
import pandas as pd
from attr import attrs,attrib

proto_filename = sys.argv[1]
standard_filename = sys.argv[2]
genes_fasta_filename = sys.argv[3]
gene_annotations_filename = sys.argv[4]
proto_out_file = sys.argv[5]
standard_out_file = sys.argv[6]

proto_df = pd.read_csv(proto_filename)
standard_df = pd.read_csv(standard_filename)

f = open(genes_fasta_filename)
gene_seqs = f.readlines()
gene_seqs = gene_seqs[1::2]
f.close()

f = open(gene_annotations_filename)
gene_annots = f.readlines()
f.close()

nt_idxs = {"A": 0, "T": 1, "C": 2, "G": 3, "X": 1}
BP_SHIFT = 6 # How far into the BP motif is the branchpoint?
THREEPRIME_SHIFT = 2 # How far into the 3'SS motif is the end of the intron?

@attrs
class IntronPWMs:
    fiveprime_pwm = attrib()
    threeprime_pwm = attrib()
    bp_pwm = attrib()
    fiveprime_cutoff = attrib()
    threeprime_cutoff = attrib()
    bp_cutoff = attrib()

@attrs
class LengthCutoffs:
    min_length = attrib()
    max_length = attrib()
    min_bp_end = attrib()
    max_bp_end = attrib()
    min_start_bp = attrib()
    max_start_bp = attrib()

@attrs
class BedLine:
    chromosome = attrib()
    start = attrib()
    end = attrib()
    name = attrib()
    strand = attrib()
    bp = attrib()

    def write(self):
        bed_tuple = (self.chromosome, str(self.start), str(self.end), self.name, "0", self.strand, str(self.bp))
        write_str = '\t'.join(bed_tuple) + '\n'
        return write_str

# Gets the PWM from a set of sequences
# Requirements: All sequences are the same length, 
# consist of ATCG, and includes at least 1 sequence.
def get_pwm(seqs):
    counts = np.array([1] * len(seqs[0]) * 4) 
    counts.resize((4, len(seqs[0])))

    for seq in seqs:
        for ii, curchar in enumerate(seq):
            counts[nt_idxs[curchar], ii] += 1

    counts = counts/counts.sum(axis=0, keepdims=True)
    return np.log(counts/0.25)

# Get sorted scores for all sequences based on the given pwm
def get_score_list(pwm, seqs):
    score_list = np.zeros(len(seqs))
    for ii, seq in enumerate(seqs): 
        align_score = 0
        for jj, curchar in enumerate(seq):
            align_score += pwm[nt_idxs[curchar], jj]
        score_list[ii] = align_score
    score_list = np.sort(score_list)
    return score_list

# Get the pwm score at a particular percentile for the list of sequences
def get_score_percentile(pwm, seqs, perc):
    score_list = get_score_list(pwm, seqs)
    return score_list[int(len(score_list) * perc)]

# Get the pwm score at each percentile for the list of sequences
def get_score_percentiles(pwm, seqs):
    score_list = get_score_list(pwm, seqs)
    percs = [0, 0.01, 0.25, 0.5, 0.75, 0.9, 0.999]
    percentile_list = []
    for perc in percs:
        percentile_list += [score_list[int(len(score_list) * perc)]]
    return percentile_list

# Get the 5'SS, BP, and 3'SS PWM's and percentile cutoffs given a 
# df that has the 5'SS, BP, and 3'SS sequences
def get_pwms(df, perc=0.1, verbose=False):
    fiveprime_seqs = list(df["5'SS"])
    bp_seqs = list(df["BP"])
    threeprime_seqs = list(df["3'SS"])

    fiveprime_pwm = get_pwm(fiveprime_seqs)
    fiveprime_cutoff = get_score_percentile(fiveprime_pwm, fiveprime_seqs, perc)
    if verbose:
        score_perc = get_score_percentiles(fiveprime_pwm, fiveprime_seqs)
        print("Getting 5'SS PWM")
        print(fiveprime_pwm)
        print(score_perc)
    bp_pwm = get_pwm(bp_seqs)
    bp_cutoff = get_score_percentile(bp_pwm, bp_seqs, perc)
    if verbose:
        score_perc = get_score_percentiles(bp_pwm, bp_seqs)
        print("Getting BP PWM")
        print(bp_pwm)
        print(score_perc)
    threeprime_pwm = get_pwm(threeprime_seqs)
    threeprime_cutoff = get_score_percentile(threeprime_pwm, threeprime_seqs, perc)
    if verbose:
        score_perc = get_score_percentiles(threeprime_pwm, threeprime_seqs)
        print("Getting 3'SS PWM")
        print(threeprime_pwm)
        print(score_perc)
    return IntronPWMs(fiveprime_pwm=fiveprime_pwm, threeprime_pwm=threeprime_pwm, bp_pwm=bp_pwm, \
        fiveprime_cutoff=fiveprime_cutoff, threeprime_cutoff=threeprime_cutoff, bp_cutoff=bp_cutoff)

def get_perc_values(values, perc):
    values = np.sort(np.array(values).astype(int))
    min_value = values[int(len(values) * perc)]
    max_value = values[int(len(values) * (1 - perc))]
    return (min_value, max_value)

# Get percentile cutoffs for lengths between start, BP, and end
def get_length_cutoffs(df, perc=0.1):
    (min_length, max_length) = get_perc_values(list(df["Length (bp)"]), perc)
    (min_bp_end, max_bp_end) = get_perc_values(list(df["BP to 3'SS (bp)"]), perc)
    (min_start_bp, max_start_bp) = get_perc_values(list(df["5'SS to BP (bp)"]), perc)
    
    return LengthCutoffs(min_length=min_length, max_length=max_length, \
        min_bp_end=min_bp_end, max_bp_end=max_bp_end, min_start_bp=min_start_bp, max_start_bp=max_start_bp)

# Assumes that cutoff > 0
# Find all indices in the sequence that match the PWM with a score at least cutoff
def find_all_matching_idxs(pwm, seq, cutoff, verbose=False):
    pwm = pwm.flatten('F')
    seq_idxs = np.array([nt_idxs[curchar] for curchar in seq])
    pwm_len = int(len(pwm)/4)
    
    # Index into the PWM is offset by this matrix
    num_repeats = int(math.ceil(len(seq)/pwm_len))
    idx_offset_full = np.tile(np.arange(pwm_len) * 4, num_repeats)
    ones_array_single = np.zeros(pwm_len)
    ones_array_single[0] = 1
    ones_array_full = np.tile(ones_array_single, num_repeats)
    idx_array = np.arange(len(seq))

    matching_idxs = np.array([])
    # Score every possible shifting of the PWM in the sequence
    for ii in range(pwm_len):
        idx_offset = idx_offset_full[:len(seq)]
        ones_array = ones_array_full[:len(seq)]
        ones_array[-(pwm_len-1):] = 0 # Make sure the overhang positions are 0. Motif matches don't count if they are not complete by the end of the sequence.

        # For each possible shift, get all PWM scores and find the idxs of the successes
        pwm_scores = pwm[seq_idxs + idx_offset]
        cur_sum = pwm_scores
        for ii in range(pwm_len - 1):
            pwm_scores = np.roll(pwm_scores, -1)
            cur_sum += pwm_scores
        cur_sum = cur_sum * ones_array # Multiply by array of 1's wherever these PWM match candidates start, 0's elsewhere
        cur_idxs = idx_array[cur_sum > cutoff]    # Get indexes where value is greater than cutoff
        matching_idxs = np.concatenate([matching_idxs, cur_idxs])
        idx_offset_full = np.roll(idx_offset_full, 1) # Rotates elements of the array 
        ones_array_full = np.roll(ones_array_full, 1)

    if verbose: # Print out the motif matches
        for idx in matching_idxs:
            print(seq[int(idx):int(idx+pwm_len)])
    
    return matching_idxs.astype(int)

# From (fiveprime_pwm, bp_pwm, threeprime_pwm), finds candidate
# 5'SS, BP, and 3'SS 
def get_introns(pwms, seq, length_cutoffs, verbose=False):
    fiveprime_candidates = find_all_matching_idxs(pwms.fiveprime_pwm, seq, pwms.fiveprime_cutoff, verbose=verbose)
    bp_candidates = find_all_matching_idxs(pwms.bp_pwm, seq, pwms.bp_cutoff, verbose=verbose)
    threeprime_candidates = find_all_matching_idxs(pwms.threeprime_pwm, seq, pwms.threeprime_cutoff, verbose=verbose)
    if verbose:
        print("Num 5'candidates:")
        print(len(fiveprime_candidates))
        print("Num BP candidates:")
        print(len(bp_candidates))
        print("Num 3'candidates:")
        print(len(threeprime_candidates))

    # These index lists are used for quickly checking whether there
    # are BP and 3'SS candidates in the right spot
    bp_idxs = np.zeros(len(seq))
    bp_idxs[bp_candidates] = 1
    threeprime_idxs = np.zeros(len(seq))
    threeprime_idxs[threeprime_candidates] = 1

    # Find intron coordinate triplets that reflect the correct spacing of the 5'SS, BP, and 3'SS. 
    intron_coords = []
    for fiveprime_idx in fiveprime_candidates:
        bp_start = fiveprime_idx + length_cutoffs.min_start_bp + BP_SHIFT
        bp_end = fiveprime_idx + length_cutoffs.max_start_bp + 1 + BP_SHIFT
        bp_end = min(bp_end, len(seq))
        if (bp_start >= bp_end) or (sum(bp_idxs[bp_start:bp_end]) == 0):
            continue
        for bp_candidate in bp_candidates:
            bp_idx = bp_candidate + BP_SHIFT
            # Check that BP is the right distance from 5'SS
            if bp_idx >= bp_start and bp_idx <= bp_end:
                threeprime_start = bp_idx + length_cutoffs.min_bp_end + THREEPRIME_SHIFT
                threeprime_end = bp_idx + length_cutoffs.max_bp_end + THREEPRIME_SHIFT
                threeprime_end = min(threeprime_end, len(seq))
                if (threeprime_start >= threeprime_end) or \
                    ((sum(threeprime_idxs[threeprime_start:threeprime_end]) == 0)):
                    continue
                # Pick the three prime idx with the shortest "feasible" distance to the BP
                closest_threeprime_idx = threeprime_end + 1
                for threeprime_candidate in threeprime_candidates:
                    threeprime_idx = threeprime_candidate + THREEPRIME_SHIFT
                    # Check that 3'SS is the right distance from the BP
                    if threeprime_idx >= threeprime_start and threeprime_idx <= threeprime_end:
                        full_start = fiveprime_idx + length_cutoffs.min_length
                        full_end = fiveprime_idx + length_cutoffs.max_length
                        full_end = min(full_end, len(seq))
                        # Check that 3'SS is the right distance from 5'SS
                        if full_start < full_end and \
                            threeprime_idx >= full_start and threeprime_idx <= full_end:
                            if threeprime_idx < closest_threeprime_idx:
                                closest_threeprime_idx = threeprime_idx
                if closest_threeprime_idx <= threeprime_end:
                    intron_coords += [(fiveprime_idx, bp_idx, closest_threeprime_idx)]
    return intron_coords

def reverse_invert(intron):
    replace_ch = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_intron = ''
    for ii in range(0, len(intron)):
        if (intron[ii] == '\n'):
            continue
        reversed_intron = replace_ch[intron[ii]] + reversed_intron
    return reversed_intron

# Assemble BED file using pre-assembled PWMs and length cutoffs from a 
# list of gene sequences and annotations
def get_candidate_beds(pwms, length_cutoffs, gene_seqs, gene_annots, one_per_gene=False, verbose=False):
    candidate_beds = []

    for ii, seq in enumerate(gene_seqs):
        if (ii % 100) == 0:
            print("Processing gene: " + str(ii) + "/" + str(len(gene_seqs)))
        seq = seq.split('\n')[0]

        # Get context for BED file coordinates, reverse sequence if needed
        gene_annot = gene_annots[ii]
        annot_items = gene_annot.split('\t')
        chromosome = annot_items[0]
        t_start = annot_items[1]
        t_end = annot_items[2]
        strand = annot_items[5].strip('\n')[0]
        if (strand == "-"):
            seq = reverse_invert(seq)
            t_start = annot_items[2]
            t_end = annot_items[1]
        t_start = int(t_start)
        t_end = int(t_end)
        gene_name = annot_items[3]

        intron_candidates = get_introns(pwms, seq, length_cutoffs)
        if verbose:
            for (five, bp, three) in intron_candidates:
                print((five, bp, three))
                print(seq[five:(five+6)])
                print(seq[(bp-BP_SHIFT):(bp-BP_SHIFT+7)])
                print(seq[(three-THREEPRIME_SHIFT):(three-THREEPRIME_SHIFT+3)])
        if verbose:
            if len(intron_candidates) > 0:
                print(len(intron_candidates))

        # To reduce the number of "fake introns", require that there is at most one per gene
        # Just choose it randomly for now
        if one_per_gene and len(intron_candidates) > 0:
            intron_candidates = [intron_candidates[int(len(intron_candidates) * random.random())]]

        # Get BED File coordinates
        cur_beds = []
        for intron_candidate in intron_candidates:
            start_loc = t_start + intron_candidate[0]
            end_loc = t_start + intron_candidate[2] + 1
            if (strand == "-"):
                start_loc = t_start - intron_candidate[2] - 1
                end_loc = t_start - intron_candidate[0]
            # Get branchpoint location relative to intron start. I checked that this gives the right branchpoint.
            branchpoint_loc = intron_candidate[1]-intron_candidate[0]-1
            cur_beds += [BedLine(chromosome=chromosome, start=start_loc, end=end_loc, \
                name=gene_name, strand=strand, bp=branchpoint_loc)] 
        candidate_beds += cur_beds

    return candidate_beds

def print_bed_file(candidate_beds, out_file):
    f = open(out_file, 'w')
    for candidate_bed in candidate_beds:
        f.write(candidate_bed.write())
    f.close()

proto_pwms = get_pwms(proto_df, perc=0.2)
length_cutoffs = get_length_cutoffs(proto_df, perc=0.2)
print(length_cutoffs)
candidate_beds = get_candidate_beds(proto_pwms, length_cutoffs, gene_seqs, gene_annots, one_per_gene=True)
print(len(candidate_beds))
print_bed_file(candidate_beds, proto_out_file)

standard_pwms = get_pwms(standard_df, perc=0.05)
length_cutoffs = get_length_cutoffs(standard_df, perc=0.05)
print(length_cutoffs)
candidate_beds = get_candidate_beds(standard_pwms, length_cutoffs, gene_seqs, gene_annots)
print(len(candidate_beds))
print_bed_file(candidate_beds, standard_out_file)

