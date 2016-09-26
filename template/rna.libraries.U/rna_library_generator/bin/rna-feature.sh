#!/bin/bash
for rna in `cat $1 | awk '{print $1}'`; do 

	grep $rna ./database/rna.txt | awk '{print $2}' > ./tmp/seq-rna.txt; 

	# secondary structure
	bash ./bin/secondary.structures.sh ./tmp/seq-rna.txt 6 $rna > ./results/check;
	
	if [[ `wc -l ./results/check | awk '{print $1}'` -lt 6 ]]; then
	
		echo "$rna" >> errors.log
	
	else
	
		cat ./results/check
		cp ./results/out ./results/$rna.ss.txt

		#VdW and polar
		bash ./bin/vdw-h.sh ./tmp/seq-rna.txt $rna > ./tmp/vdw.txt
		awk '{a[NR]=$2; b[NR]=$3} END{printf "%s\t", "'$rna'"; for(i=1;i<=NR;i++){printf "%4.2f\t",a[i]} 
		printf "\n%s\t", "'$rna'"; for(i=1;i<=NR;i++){printf "%4.2f\t",b[i]} printf "\n";}' ./tmp/vdw.txt

		bash ./bin/vdw-h.2.sh ./tmp/seq-rna.txt $rna > ./tmp/vdw.txt
		awk '{a[NR]=$2; b[NR]=$3} END{printf "%s\t", "'$rna'"; for(i=1;i<=NR;i++){printf "%4.2f\t",a[i]} 
		printf "\n%s\t", "'$rna'"; for(i=1;i<=NR;i++){printf "%4.2f\t",b[i]} printf "\n";}' ./tmp/vdw.txt 

	fi

done
