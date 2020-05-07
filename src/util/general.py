from config import DATABASE_PATH
import pandas as pd

def reverse_invert(seq):
	replace_ch = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	reversed_seq = ''
	for ii in range(0, len(seq)):
		if (seq[ii] == '\n'):
			continue
		reversed_seq = replace_ch[seq[ii]] + reversed_seq
	return reversed_seq

def get_rpkms():
	gene_exp_file = DATABASE_PATH + 'gene_expression_Hickman.csv'
	gene_exp = pd.read_csv(gene_exp_file)
	genes = gene_exp['Gene']
	rpkms_list = gene_exp['average.rpkm']
	rpkms = {}
	for ii, gene in enumerate(genes):
		rpkms[gene] = rpkms_list[ii]
	return rpkms