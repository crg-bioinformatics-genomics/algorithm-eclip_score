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
#!/bin/bash

for rna in `cat $1 | awk '{print $1}'`; do 

grep $rna ./database/rna.txt | awk '{print $2}' > ./tmp/seq-rna.txt; 

# secondary structure
sh secondary.structures.sh ./tmp/seq-rna.txt 6 $rna; 
cp ./results/out ./results/$rna.ss.txt

#VdW and polar
sh vdw-h.sh ./tmp/seq-rna.txt $rna > ./tmp/vdw.txt
awk '{a[NR]=$2; b[NR]=$3} END{printf "%s\t", "'$rna'"; for(i=1;i<=NR;i++){printf "%4.2f\t",a[i]} 
printf "\n%s\t", "'$rna'"; for(i=1;i<=NR;i++){printf "%4.2f\t",b[i]} printf "\n";}' ./tmp/vdw.txt

sh vdw-h.2.sh ./tmp/seq-rna.txt $rna > ./tmp/vdw.txt
awk '{a[NR]=$2; b[NR]=$3} END{printf "%s\t", "'$rna'"; for(i=1;i<=NR;i++){printf "%4.2f\t",a[i]} 
printf "\n%s\t", "'$rna'"; for(i=1;i<=NR;i++){printf "%4.2f\t",b[i]} printf "\n";}' ./tmp/vdw.txt 
done
