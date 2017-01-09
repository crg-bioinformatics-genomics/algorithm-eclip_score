#!/bin/bash
# bash omixcore.sh <random_number>  <email_address> <title> <type> <protein file><mode> <rnafolder>
echo $1 $2 $3 $4 $5 $6 $7

set -e
#set -o pipeline
mode=$6
rnafolder=$7
script_folder=$(pwd)
if [ -s tmp/$1 ]; then
	rm -fr tmp/$1
fi
cp -r template tmp/$1
cd   tmp/$1
cp $5 data/inseq.fasta
signature_prediction=`python RBPprofiler.py`
rbp=$(echo $signature_prediction | awk '{print $NF}')

if [[ $rbp > 0.5 ]] || [[ $rbp = 0.5 ]]
then
	# modifies fasta files into 2 columns
	# memo: outfile is the protein file
	./flip -u $5
	sed 's/[\| | \\ | \/ | \* | \$ | \# | \? | \! ]/./g' $5 | awk '(length($1)>=1)' | awk '($1~/>/){gsub(" ", "."); printf "\n%s\t", $1} ($1!~/>/){printf "%s", toupper($1)}' | awk '(NF==2)' | head -1 | sed 's/>\.//g;s/>//g' | awk '{print substr($1,1,12)"_"NR, toupper($2)}' >  ./protein/outfile

	if [[ ! -s "./protein/outfile" ]]; then
	        echo "Protein file has not been created."
		exit 1
	fi




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
			mv outfile.fr.prot.lib ../../outputs/protein.lib
			rm outfile.fr.txt
		cd ../
	cd ../

	if [[ $mode == "custom" ]]
	then
		echo "Rna sequence processing and fragmentation"
		if [ ! -s rna/rna_seqs_oneline/ ]; then
			mkdir rna/rna_seqs_oneline/
		fi
		if [ ! -s rna/rna_seqs_fragmented ]; then
			mkdir rna/rna_seqs_fragmented
		fi
		for i in `ls $rnafolder`
		do
			sed 's/[\| | \\ | \/ | \* | \$ | \# | \? | \! ]/./g' $rnafolder/$i | awk '(length($1)>=1)' | awk '($1~/>/){gsub(" ", "."); printf "\n%s\t", $1} ($1!~/>/){gsub(/[Uu]/, "T", $1); printf "%s",toupper($1)}' | awk '(NF==2)' | head -1 | sed 's/>\.//g;s/>//g' | awk '{print substr($1,1,12)"_"NR, $2}' >  ./rna/rna_seqs_oneline/$i
		done
		cd fragmentation
			echo "Custom Rna fragmentation"
			for i in `ls ../rna/rna_seqs_oneline/`
				do
					bash rna.job.cutter.sh ../rna/rna_seqs_oneline/$i > ../rna/rna_seqs_fragmented/$i.rna.fragm.seq;
				done
		cd ..

		cd ./rna/rna_seqs_fragmented

			if [ ! -s ../../rna.libraries.U/rna_seqs ]; then
				mkdir ../../rna.libraries.U/rna_seqs
			fi
			n=0
			for i in *
			do
  			if [ $((n+=1)) -gt 10 ]; then
    			n=1
  			fi
  			todir=../../rna.libraries.U/rna_seqs/$n
  			[ -d "$todir" ] || mkdir "$todir"
  			mv "$i" "$todir"
		done
		cd ../../

		echo "Custom rna library generation"

		cd rna.libraries.U/
			if [ ! -s outs ]; then
				mkdir outs
			fi
			bash superjob.sh
		cd ..
		rna_lib_folder=$(pwd | awk '{print $0"/rna.libraries.U/outs/"}')
		python library_checker.py $rna_lib_folder
		cp -r $rna_lib_folder outputs
		mv outputs/outs/ outputs/rna_libs

	fi
	if [[ $mode == "lincrnas" ]]
	then
		rna_lib_folder=$(echo $script_folder | awk '{print $0"/lincrnas/"}')
	fi
	# PROTEIN + RNA INTERACTIONS #################################################################################

	echo "protein and RNA interaction computing"
	date +"%m-%d-%y %r"

	cd interactions.U/
		cp ../outputs/protein.lib combine_parallel/prot/prot.lib

			cd combine_parallel/
				if [ ! -s pre-compiled ]; then
					mkdir pre-compiled/
				fi

				#all out.merged are saved to pre-compiled but then is distributed to folders in order to be paralelised
				python multiplier.py 10 "prot" $rna_lib_folder
				#echo "# protein / rna / raw score / dp " > ../../outputs/interactions.$1.$3.txt
				#cat pre-compiled/* >> ../../outputs/interactions.$1.$3.txt
				#cat pre-compiled/* > pre-compiled/out.merged.txt
			cd ../../

		cd filter

				echo "Global Score computing"
				date +"%m-%d-%y %r"
				NUMJOBS=10
				n=0
				for i in `ls ../interactions.U/combine_parallel/pre-compiled/out.merged.*`
				do
				  if [ $((n+=1)) -gt $NUMJOBS ]; then
				    n=1
				  fi
				  todir=../interactions.U/combine_parallel/pre-compiled/$n
				  [ -d "$todir" ] || mkdir "$todir"
				  mv "$i" "$todir"

 				done

				for j in `ls ../interactions.U/combine_parallel/pre-compiled/`
				do
					( cp -r template dir.$j
					for i in `ls ../interactions.U/combine_parallel/pre-compiled/$j`
					do
						cd dir.$j
						prot_rna=$(echo $i | awk -F 'out.merged.' '{print $2}' |  awk -F '.txt' '{print $1}')
						rna=$(echo $prot_rna | awk -F '-' '{print $2}')

						awk '{print "protein_"$1,"rna_"$2,$3,$4}' ../../interactions.U/combine_parallel/pre-compiled/$j/$i > interactions.$prot_rna.txt
						bash gs_network.sh interactions.$prot_rna.txt $rna

						paste tmp/names.txt tmp/output.txt | awk '{print $1, $NF}'> $prot_rna.processed.txt
						#bash start.sh interactions.$prot_rna.txt > $prot_rna.processed.txt
						echo $rna $(awk '{printf "%.2f\n", ($2+1)/2}' $prot_rna.processed.txt) >> filter.processed.txt
						cd ..
					done ) &
				done
				wait
				cat dir.*/filter.processed.txt > ../outputs/filter.processed.txt
				rm -fr dir.*
			cd ..
			rm -fr interactions.U/combine_parallel/pre-compiled/*
			rm -fr images
			rm -fr data



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
