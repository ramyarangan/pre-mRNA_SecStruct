import sys
import numpy as np 
import pandas as pd
from attr import attrs,attrib

standard_filename = sys.argv[1]
score_filename = sys.argv[2]

standard_df = pd.read_csv(standard_filename)

nt_idxs = {"A": 0, "T": 1, "C": 2, "G": 3, "X": 1}
BP_SHIFT = 6 # How far into the BP motif is the branchpoint?
THREEPRIME_SHIFT = 2 # How far into the 3'SS motif is the end of the intron?

@attrs
class IntronPWMs:
    fiveprime_pwm = attrib()
    threeprime_pwm = attrib()
    bp_pwm = attrib()

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

def get_score(pwm, seq):
    align_score = 0
    for jj, curchar in enumerate(seq):
        align_score += pwm[nt_idxs[curchar], jj]
    return align_score

# Get the 5'SS, BP, and 3'SS PWM's and percentile cutoffs given a 
# df that has the 5'SS, BP, and 3'SS sequences
def get_pwms(df, perc=0.1, verbose=False):
    fiveprime_seqs = list(df["5'SS"])
    bp_seqs = list(df["BP"])
    threeprime_seqs = list(df["3'SS"])

    fiveprime_pwm = get_pwm(fiveprime_seqs)
    bp_pwm = get_pwm(bp_seqs)
    threeprime_pwm = get_pwm(threeprime_seqs)

    return IntronPWMs(fiveprime_pwm=fiveprime_pwm, threeprime_pwm=threeprime_pwm, bp_pwm=bp_pwm)

def write_length_pwm_scores(df):
    pwms = get_pwms(df)
    chrs = list(df["Chromosome"])
    starts = list(df["Start"])
    ends = list(df["End"])
    lengths = np.array(list(df["Length (bp)"])).astype(int)
    start_lengths = np.array(list(df["5'SS to BP (bp)"])).astype(int)
    end_lengths = np.array(list(df["BP to 3'SS (bp)"])).astype(int)
    fiveprime_seqs = list(df["5'SS"])
    bp_seqs = list(df["BP"])
    threeprime_seqs = list(df["3'SS"])
    seqs = list(df["sequence"])

    f = open(score_filename, 'w')
    for ii, seq in enumerate(seqs):
        fiveprime_score = get_score(pwms.fiveprime_pwm, fiveprime_seqs[ii])
        bp_score = get_score(pwms.bp_pwm, bp_seqs[ii])
        threeprime_score = get_score(pwms.threeprime_pwm, threeprime_seqs[ii])

        full_len = lengths[ii]
        start_len = start_lengths[ii]
        end_len = end_lengths[ii]

        chr_tag = chrs[ii] + ":" + str(starts[ii]) + "-" + str(ends[ii])

        f.write("%s %f %f %f %d %d %d\n" % (chr_tag, fiveprime_score, bp_score, \
            threeprime_score, full_len, start_len, end_len))

    f.close()


write_length_pwm_scores(standard_df)
