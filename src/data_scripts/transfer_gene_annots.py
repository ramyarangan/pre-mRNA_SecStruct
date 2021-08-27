import sys

input_datfile = sys.argv[1]
input_annotfile = sys.argv[2]
output_datfile = sys.argv[3]
offset_5ss = sys.argv[4]
offset_3ss = sys.argv[5]

f = open(input_datfile)
dat_lines = f.readlines()
f.close()

f = open(input_annotfile)
annot_lines = f.readlines()
f.close()

f = open(output_datfile, 'w')

chr_dict = {}
for ii in range(int(len(annot_lines)/2)):
	annot_line = annot_lines[2 * ii]
	annot_items = annot_line.strip('\n').split(' ')[1].split('\t')
	if annot_items[0] not in chr_dict.keys():
		chr_dict[annot_items[0]] = []
	cur_annot_items = annot_items[1:]
	cur_annot_items[0] = int(cur_annot_items[0])
	cur_annot_items[1] = int(cur_annot_items[1])
	chr_dict[annot_items[0]] += [cur_annot_items]

for ii in range(int(len(dat_lines)/2)):
	dat_line = dat_lines[2 * ii].strip('\n')
	dat_bp = dat_line.split(' ')[0]
	dat_items = dat_line.split(' ')[1].split('\t')
	annot_item = chr_dict[dat_items[0]]
	intron_start = int(dat_items[1])
	intron_end = int(dat_items[2])
	found_annot = None
	for chr_annot in annot_item:
		if chr_annot[0] < intron_start and chr_annot[1] > intron_end:
			found_annot = chr_annot
			break
		if chr_annot[0] >= intron_start and chr_annot[1] <= intron_end:
			found_annot = chr_annot
			break
	if found_annot is None:
		print("FAIL")
		continue
	found_annot += [offset_5ss]
	found_annot += [offset_3ss]
	found_annot = [str(x) for x in found_annot]
	f.write('%s %s\t%s\n' % (dat_bp, dat_items[0], '\t'.join(found_annot)))
	f.write('%s\n' % dat_lines[2 * ii + 1].strip('\n'))

f.close()