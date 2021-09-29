"""
Writes all secondary structure predictions to a file

Example usage: python analyses/write_secstructs.py standard_allsize util/tmp/intron_mfe.txt
"""

import os
import argparse
from config import PATH_VARNA_JAR
from util.features_db import * 

parser = argparse.ArgumentParser(description='Parameters for writing secondary structures')
parser.add_argument('intron_class', type=str, help='Intron class')
parser.add_argument('secstruct_file', type=str, help='Output file')
parser.add_argument('--package', type=str, help='secondary structure package', default='contrafold')

args = parser.parse_args()
intron_class = args.intron_class
secstruct_file = args.secstruct_file
package = args.package

secstruct_options = {'secstruct_pkg': package, 
					'secstruct_type': 'mfe', 
					'verbose': False,
					'force_eval': False
					}

check_update_secstruct(intron_class, secstruct_options)
all_introns = build_intron_set(intron_class, secstruct_options)

f = open(secstruct_file, 'w')

for intron in all_introns.introns:
	f.write("%s\n" % intron.print_string())
	f.write("%s\n" % intron.seq)
	f.write("%s\n" % intron.mfe)

f.close()