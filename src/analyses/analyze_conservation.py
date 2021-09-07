import argparse


for zipper_stem in zipper_stem_file_data:
	[seq1, seq2, secstruct, dG] = zipper_stem


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Parameters for processing intron alignment data')
	parser.add_argument('alignment_dir', type=str, help='Path to directory storing alignments in Stockholm format')
	parser.add_argument('intron_class', type=str, help='Path to directory storing alignments in Stockholm format')
	parser.add_argument('--do_zipper', default=False, action='store_true', \
		 help='Do stats on zipper stems')	
	args = parser.parse_args()

	alignment_dir = args.alignment_dir
	intron_class = args.intron_class
	do_zipper = args.do_zipper

	get_intron_aln(intron_seq, alignment_dir)