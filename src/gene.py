from config import DATABASE_PATH
from util.general import reverse_invert

STOP_CODONS = ['TAG', 'TAA', 'TGA']

@attrs
class Gene:
	chr_num = attrib()
	chr_pos = attrib()
	refseq_name = attrib()
	strand_dir = attrib()
	seq = attrib()

class GeneSet:
	def __init__(self, prune_ORFs=True):
		gene_seq_file = DATABASE_PATH + 'genes/all_genes.fasta'
		gene_annot_file = DATABASE_PATH + 'genes/gene_annotations.bed'

		f = open(gene_seq_file)
		gene_seqs = f.readlines()
		f.close()

		f = open(gene_annot_file)
		gene_annots = f.readlines()
		f.close()

		genes = []
		for ii in range(len(gene_annots)):
			gene_annot = gene_annots[ii][:-1]
			seq = gene_seqs[ii * 2 + 1][:-1]
			info_items = gene_annot.split('\t')
			chr_num = info_items[0]
			chr_pos = (int(info_items[1]), int(info_items[2]))
			refseq_name = info_items[3]
			strand_dir = info_items[5]
			if strand_dir == '-':
				seq = reverse_invert(seq)
			gene = Gene(chr_num=chr_num, chr_pos=chr_pos, refseq_name=refseq_name, 
						strand_dir=strand_dir, seq=seq)
			genes += [gene]
		if prune_ORFs:
			genes = prune_ORFs(genes)
		self.genes = genes 

	def get_genes_dict(self):
		genes_dict = {} # From refseq name to gene object
		for gene in self.genes:
			genes_dict[gene.refseq_name] = gene
		return genes_dict

	def prune_ORFs(genes):
		# Clean genes list.. keep only sequences that are ORFs
		new_genes = []
		for gene in genes:
			if gene.seq[:3] != 'ATG':
				continue
			if gene.seq[-3:] not in STOP_CODONS:
				continue
			new_genes += [gene]
		return new_genes

