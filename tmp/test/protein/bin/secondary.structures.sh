#***************************************************************************
#                             catRAPID 
#                             -------------------
#    begin                : Nov 2010
#    copyright            : (C) 2010 by G.G.Tartaglia
#    author               : : Tartaglia
#    date                 : :2010-11-25
#    id                   : caRAPID protein-RNA interactions
#    email                : gian.tartaglia@crg.es
#**************************************************************************/
#/***************************************************************************
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU General Public License as published by  *
# *   the Free Software Foundation; either version 2 of the License, or     *
# *   (at your option) any later version.                                   *
# *                                                                         *
# ***************************************************************************/
cp ./$1 ./tmp/seq.txt

# creates $2 models

RNAsubopt -e 1 -s < ./tmp/seq.txt | head -10 | sort -n -k 2 | awk '{print $1}' > ./results/out

#
for((ix=2;ix<=$2+1;ix++))
do
awk '((NR==1)||NR=='$ix')' ./results/out > ./results/3.out

# generates the coordinates
RNAplot -o xrna < ./results/3.out
# cleans the file
awk '(NR>4)' rna.ss   > ./results/clean.txt;

# creates the contact profile
awk '{a[NR]=$3; b[NR]=$4;} END{for(i=1;i<=NR;i++){k=0; for(j=1;j<=NR;j++){ if((i!=j)&&(sqrt((a[i]-a[j])^2+(b[i]-b[j])^2))<50){k++;}}  print i,k/50}}' ./results/clean.txt > ./results/4.out
awk 'BEGIN{printf "%s\t", "'$3'"; } {printf "%4.2f\t", $2} END{printf "\n"}' ./results/4.out;
done; 

