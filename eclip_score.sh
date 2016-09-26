#!/bin/bash
# bash runcatrapid.fragments.sh <random_number> <email_address> <title> <type> <sizeRNAs> <protein file> <RNA file> <extras>
echo $1 $2 $3 $4 $5 $6 $7 $8

set -e
#set -o pipeline

if [ -s tmp/$1 ]; then
	rm -fr tmp/$1
fi
cp -r template tmp/$1
cd   tmp/$1

# modifies fasta files into 2 columns
# memo: outfile is the protein file
./flip -u $6
sed 's/[\| | \\ | \/ | \* | \$ | \# | \? | \! ]/./g' $6 | awk '(length($1)>=1)' | awk '($1~/>/){gsub(" ", "."); printf "\n%s\t", $1} ($1!~/>/){printf "%s", toupper($1)}' | awk '(NF==2)' | head -1 | sed 's/>\.//g;s/>//g' | awk '{print substr($1,1,12)"_"NR, toupper($2)}' >  ./protein/outfile
./flip -u $7
sed 's/[\| | \\ | \/ | \* | \$ | \# | \? | \! ]/./g' $7 | awk '(length($1)>=1)' | awk '($1~/>/){gsub(" ", "."); printf "\n%s\t", $1} ($1!~/>/){gsub(/[Uu]/, "T", $1); printf "%s",toupper($1)}' | awk '(NF==2)' | head -1 | sed 's/>\.//g;s/>//g' | awk '{print substr($1,1,12)"_"NR, $2}' >  ./rna/outfile2

if [[ ! -s "./protein/outfile" || ! -s "./rna/outfile2" ]]; then
        echo "A file has not been created."
	exit 1
fi

echo "# Reference Sequences" >  ./outputs/yoursequences.$1.$3.txt
cat ./protein/outfile  	    >>  ./outputs/yoursequences.$1.$3.txt
cat ./rna/outfile2	        >>  ./outputs/yoursequences.$1.$3.txt

num1=`cat ./rna/outfile2     | awk '{print NR}'`
num2=`cat ./protein/outfile  | awk '{print NR}'`

echo "sequence processing and fragmentation"
date +"%m-%d-%y %r"

cd ./rna

	echo "# SERVER RUNNING...PLEASE REMEMBER TO REFRESH THIS PAGE FROM TIME TO TIME..." >  ../outputs/library.$1.$3.txt;
	nt=`cat outfile2           | awk '(NR==1){s=length($2)/50; if(s<=25){s=25} if(s>=750){s=750} print int(s)}'`
	bash ../fragmentation/cutter.sh outfile2 $nt >  outfile2.fr
cd ../

cd protein
	pr=`cat ../protein/outfile | awk '(NR==1){s=length($2)/50; if(s<=25){s=25} if(s>=375){s=375} print int(s)}'`
	bash ../fragmentation/cutter.sh outfile  $pr >  outfile.fr
cd ../

	# RNA LIBRARY GENERATION ######################################################################################
echo "protein and RNA libraries generation"
date +"%m-%d-%y %r"
	# run cases one by one
cd rna.libraries.U
	cd rna_library_generator/
		sed 's/|/_/g' ../../rna/outfile2.fr > outfile2.fr.txt
		bash rungenerator.rna.sh outfile2.fr.txt
		mv outfile2.fr.rna.lib ../../outputs/rna.library.$1.$3.txt
		rm outfile2.fr.txt
	cd ../
cd ../

wc rna/outfile2.fr  | awk '($1==0){print "Sorry! no fragments found!"}' >> outputs/library.$1.$3.txt

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
	cp ../outputs/rna.library.$1.$3.txt combine_parallel/rna/rna.lib

		cd combine_parallel/
			mkdir pre-compiled/
			bash parallelize.sh 8 prot/prot.lib rna/rna.lib
			echo "# protein / rna / raw score / dp " > ../../outputs/interactions.$1.$3.txt
			cat pre-compiled/out.merged.txt >> ../../outputs/interactions.$1.$3.txt
			rm prot/prot.lib rna/rna.lib

		cd ..
cd ..


wc protein/outfile.fr  | awk '($1==0){print "Sorry! no fragments found!"}' >> outputs/interactions.$1.$3.txt




########################
####### FILTER #########
########################


##### filter for fragments

## N.B. With UNIFORM option, IDs are removed and just fragment coordinates are shown;
## For this reason, we need to add "protein_" and "rna_" before the coordinates.

cd filter

if [[ "$4" == "uniform" ]]; then

	echo "Global Score computing"
	date +"%m-%d-%y %r"

	awk '{print "protein_"$1,"rna_"$2,$3,$4}' ../interactions.U/combine_parallel/pre-compiled/out.merged.txt > interactions.$1.$3.txt
	bash start.sh interactions.$1.$3.txt -1 > processed.txt
	awk '{printf "%.2f\n", ($2+1)/2}' processed.txt > ../outputs/filter.tmp
fi

cd ..


######################################

# MAKES STATISTICS   ########################################################################################
#

echo "statistics and table"
date +"%m-%d-%y %r"

if [[ "$4" == "weighted" ]] || [[ "$8" != "signalLoc" ]]; then
cd ./statistics

	if [ "$opt" = "1" ] ; then # uniform
		bash bin/run.2.sh $1 $3
	fi

	if [ "$opt" = "0" ] ; then # weighted
		bash bin/run.4.sh $1 $3
	fi
	# rm cases
	# rm *pdf

cd ..
awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td>%s</td>\n\t\t<td>%.2f</td>\n\t\t<td>%d</td>\n\t\t<td>%.2f</td>\n\t</tr>\n", $1, $2, $3, $4, $5*100, $6}END{print """</tbody>"""}' outputs/best.prot.rna.txt > outputs/statistics.$1.$3.html
fi

if [[ "$4" == "weighted" ]]; then
php ../../index.Z.html > outputs/index.Z.tmp.html
fi
if [[ "$4" == "uniform" ]]; then
php ../../index.mat.html > outputs/index.mat.tmp.html
fi

echo "spaghetti done!"
date +"%m-%d-%y %r"

#awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td>%s</td>\n\t\t<td>%.2f</td>\n\t\t<td>%d</td>\n\t\t<td>%.2f</td>\n\t</tr>\n", $1, $2, $3, $4, $5*100, $6}END{print """</tbody>"""}' outputs/best.prot.rna.txt > outputs/statistics.$1.$3.html
