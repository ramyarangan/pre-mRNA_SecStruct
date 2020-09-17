#!/bin/bash

seq_file=$1
out_pref=$2
total_introns=$(wc -l $seq_file | awk '{ print $1 }')
echo $total_introns

for i in $(seq 1 2 $total_introns);
do
	echo $i
	line1=$(head -n $(($i+1)) $seq_file | tail -n 2 | head -n 1)
	line2=$(head -n $(($i+1)) $seq_file | tail -n 1)
	echo ">"$line1 >> tmp.fasta
	echo $line2 >> tmp.fasta

	# RNAstructure
	#echo ">"$line1 >> $out_pref"RNAstructure_mfe.dat"
	#Fold tmp.fasta - -mfe -k | tail -n 1 >> $out_pref"RNAstructure_mfe.dat"

	# echo ">"$line1 >> $out_pref"RNAstructure_ens.dat"
	#stochastic tmp.fasta test.ct -e 1000 --sequence > tmp.out
	#ct2dot test.ct -1 - | awk 'NR%2==1' >> $out_pref"RNAstructure_ens.dat"
	#rm test.ct
	#rm tmp.out

	# Vienna
	echo ">"$line1 >> $out_pref"Vienna_mfe.dat"
	RNAfold -i tmp.fasta | tail -n 1 | awk '{ print $1 }' >> $out_pref"Vienna_mfe.dat"

	echo ">"$line1 >> $out_pref"Vienna_ens.dat"
	RNAsubopt --o=test.dat --stochBT_en=1000 < tmp.fasta
	cat test.dat | tail -n +3 | awk '{ print $1 }' >> $out_pref"Vienna_ens.dat"
	rm test.dat
	
	# Contrafold
	#echo ">"$line1 >> $out_pref"contrafold_mfe.dat"
	#contrafold predict tmp.fasta | tail -n 1 >> $out_pref"contrafold_mfe.dat"

	rm tmp.fasta
	rm *.ps
done
