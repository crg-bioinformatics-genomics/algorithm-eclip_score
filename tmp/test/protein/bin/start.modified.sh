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
# please read and install tester.cpp

#!/bin/bash
set -e

# takes the sequences (label - sequence)
echo  $1 $2 > ../database/protein.txt

# r=21; 

awk '{print "'$1'", $1}' ../rna.list.txt | sed 's/ENSU/ENST/g' | sed 's/>//g' > ../database/interactions.rna.50.750.txt 
# $1 preceeds $3
#echo $1 $3 > ./database/interactions.rna.txt



# computes 10 features (concatenated) secondary + polarity + hydro for rna and proteins
bash bin/protein-feature.sh    ../database/protein.txt   > ../database/protein.dat

# normalizes the lengths
bash bin/adaptator.sh        ../database/protein.dat 751      > ../database/protein.750.dat

# computes fouriers' coefficients
bash bin/fourier.line.sh     ../database/protein.750.dat      > ../database/fourier.protein.750.dat

# creates the inputs
# bait is the protein
ln=`wc ../rna.list.txt | awk '{print $1}'`
awk '{for(i=2;i<=51;i++){printf "\t%4.2f", $i} printf "\n"}'   ../database/fourier.protein.750.dat | awk '{a[NR]=$0} END{for(i=1;i<='$ln';i++){for(j=1;j<=10;j++){print a[j]} }}' > ../data/bait.posi

cp  ../data/*     ./data/
cp  ../database/* ./database/

# computes the score
./bin/tester 1 1 50  >    ../tmp/log.tmp

# calculates the score
bash bin/processer.sh positives.txt > ../tmp/score.txt

paste ../database/interactions.rna.50.750.txt ../tmp/score.txt | sed 's/ENSU/ENST/g'

