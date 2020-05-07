def reverse_invert(seq):
    replace_ch = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_seq = ''
    for ii in range(0, len(seq)):
        if (seq[ii] == '\n'):
            continue
        reversed_seq = replace_ch[seq[ii]] + reversed_seq
    return reversed_seq