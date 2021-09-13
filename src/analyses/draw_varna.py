import os
from arnie.mfe import mfe
import argparse
from config import PATH_VARNA_JAR
from util.features_db import * 

def get_varna_command(sequence, secstruct, mapping_data):
	varna_cmd = 'java -cp %s  fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN "%s" -colorMapMin "-0.001" -colorMapMax "2.0" -colorMapStyle "%s" -colorMap "%s"' % \
		(PATH_VARNA_JAR, sequence, secstruct, "-0.001:#C0C0C0,-0.0005:#FFFFFF;0.1:#FFFFFF,0.8:#FF8800;1:#FF0000", ','.join([str(x) for x in mapping_data]))
	return varna_cmd

def draw_varna(sequence, secstruct, shape_signal=[]):
	if len(shape_signal) == 0:
		shape_signal = [0] * len(sequence)
	varna_command = get_varna_command(sequence, secstruct, shape_signal)
	os.system(varna_command)



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Parameters for drawing VARNA figure for intron')
	parser.add_argument('ensembl_name', type=str, help='Ensembl name to plot')
	args = parser.parse_args()

	ensembl_name = args.ensembl_name

	secstruct_options = {'secstruct_pkg': 'Vienna', 
						'secstruct_type': 'mfe', 
						'verbose': False,
						'force_eval': False
						}
	intron_class = 'standard_allsize'
	check_update_secstruct(intron_class, secstruct_options)
	all_introns = build_intron_set(intron_class, secstruct_options)
	intron = all_introns.get_intron_by_ensembl_name(ensembl_name)
	draw_varna(intron.seq, intron.mfe)