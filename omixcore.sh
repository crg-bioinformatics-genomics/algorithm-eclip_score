#!/bin/bash
# bash omixcore.sh <random_number> <email_address> <title> <type> <protein file>
echo $1 $2 $3 $4 $5

set -e
#set -o pipeline

script_folder=$(pwd)
if [ -s tmp/$1 ]; then
	rm -fr tmp/$1
fi
cp -r template tmp/$1
cd   tmp/$1

signature_prediction=`python RBPprofiler.py`
rbp=$(echo $signature_prediction | awk '{print $NF}')

if (( rbp > =0.5))
then
	# modifies fasta files into 2 columns
	# memo: outfile is the protein file
	./flip -u $5
	sed 's/[\| | \\ | \/ | \* | \$ | \# | \? | \! ]/./g' $5 | awk '(length($1)>=1)' | awk '($1~/>/){gsub(" ", "."); printf "\n%s\t", $1} ($1!~/>/){printf "%s", toupper($1)}' | awk '(NF==2)' | head -1 | sed 's/>\.//g;s/>//g' | awk '{print substr($1,1,12)"_"NR, toupper($2)}' >  ./protein/outfile

	if [[ ! -s "./protein/outfile" ]]; then
	        echo "Protein file has not been created."
		exit 1
	fi

	echo "# Reference Sequences" >  ./outputs/yoursequences.$1.$3.txt
	cat ./protein/outfile  	    >>  ./outputs/yoursequences.$1.$3.txt


	num2=`cat ./protein/outfile  | awk '{print NR}'`

	echo "Protein sequence processing and fragmentation"
	date +"%m-%d-%y %r"

	cp protein/outfile fragmentation/outfile
	cd fragmentation
		bash protein.job.cutter.sh outfile > ../protein/outfile.fr
	cd ../

		# RNA LIBRARY GENERATION ######################################################################################
	echo "Protein library generation"
	date +"%m-%d-%y %r"

	cd protein.libraries.U
		cd protein_library_generator/
			sed 's/|/_/g' ../../protein/outfile.fr > outfile.fr.txt
			bash rungenerator.protein.sh outfile.fr.txt
			mv outfile.fr.prot.lib ../../outputs/protein.library.$1.$3.txt
			rm outfile.fr.txt
		cd ../
	cd ../


	# PROTEIN + RNA INTERACTIONS #################################################################################

	echo "protein and RNA interaction computing"
	date +"%m-%d-%y %r"

	cd interactions.U/
		cp ../outputs/protein.library.$1.$3.txt combine_parallel/prot/prot.lib

			cd combine_parallel/
				if [ ! -s pre-compiled ]; then
					mkdir pre-compiled/
				fi

				rna_lib_folder=$(echo $script_folder | awk '{print $0"/lincrnas/"}')
				python multiplier.py 10 "prot" $rna_lib_folder
				echo "# protein / rna / raw score / dp " > ../../outputs/interactions.$1.$3.txt
				cat pre-compiled/* >> ../../outputs/interactions.$1.$3.txt
				#cat pre-compiled/* > pre-compiled/out.merged.txt
			cd ../../

		cd filter

				echo "Global Score computing"
				date +"%m-%d-%y %r"

				for i in `ls ../interactions.U/combine_parallel/pre-compiled/`; do
					prot_rna=$(echo $i | awk -F '.' '{print $3}')
					rna=$(echo $prot_rna | awk -F '-' '{print $2}')
					awk '{print "protein_"$1,"rna_"$2,$3,$4}' ../interactions.U/combine_parallel/pre-compiled/$i > interactions.$prot_rna.txt
					bash start.sh interactions.$prot_rna.txt > $prot_rna.processed.txt
					echo $rna $(awk '{printf "%.2f\n", ($2+1)/2}' $prot_rna.processed.txt) >> ../outputs/filter.processed.txt
				done


			cd ..




	wc protein/outfile.fr  | awk '($1==0){print "Sorry! no fragments found!"}' >> outputs/interactions.$1.$3.txt







	######################################

	# MAKES STATISTICS   ########################################################################################
	#

	echo "statistics and table"
	date +"%m-%d-%y %r"


	if [[ "$4" == "uniform" ]]; then
	php ../../index.mat.html > outputs/index.mat.tmp.html
	fi

	echo "spaghetti done!"
	date +"%m-%d-%y %r"

	#awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td>%s</td>\n\t\t<td>%.2f</td>\n\t\t<td>%d</td>\n\t\t<td>%.2f</td>\n\t</tr>\n", $1, $2, $3, $4, $5*100, $6}END{print """</tbody>"""}' outputs/best.prot.rna.txt > outputs/statistics.$1.$3.html
else
	php ../../index.nrbp.html > outputs/index.mat.tmp.html
fi
