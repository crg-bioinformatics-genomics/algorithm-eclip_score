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

# creates the datasets
echo "remove" >   ./data/bait.$5
cp ./data/bait.$5 ./data/prey.$5
rm ./data/bait.$5 ./data/prey.$5
cat $2 $3 > ./tmp/data.dat
cat $1 | awk '{print $1"+"$2}' > ./tmp/plus.dat
#
for j in `more ./tmp/plus.dat `; do
#finds
A=`echo $j | sed 's/+/ /g' | awk '{print $1}'`;  
B=`echo $j | sed 's/+/ /g' | awk '{print $2}'`; 
grep -w $A ./tmp/data.dat | awk '{for(i=2;i<='$4'+1;i++){printf "\t%4.2f", $i} printf "\n"}' > ./tmp/tmp.1
grep -w $B ./tmp/data.dat | awk '{for(i=2;i<='$4'+1;i++){printf "\t%4.2f", $i} printf "\n"}' > ./tmp/tmp.2 
#double-checks
AA=`wc ./tmp/tmp.1 | awk '{print $1}'`
BB=`wc ./tmp/tmp.2 | awk '{print $1}'`
awk '('$AA'==10)&&('$BB'==10)&&(NF>1)' ./tmp/tmp.1  >> ./data/bait.$5
awk '('$AA'==10)&&('$BB'==10)&&(NF>1)' ./tmp/tmp.2  >> ./data/prey.$5
done
