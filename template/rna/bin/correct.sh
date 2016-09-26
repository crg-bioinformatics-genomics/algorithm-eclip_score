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
#bin_path=`which start.modified.sh`
#bin_path=`dirname $bin_path`
params_path=./params



# calculates the score
cat          ./tmp/score.txt      > ./tmp/analysis.txt
V=`cat tmp/analysis.txt`
cat $params_path/positives.txt    >> ./tmp/analysis.txt
P=`awk '(NR==1){s=$1} (NR>1){if($1>s){k++}} END{print k/(NR-1)}' tmp/analysis.txt`
#echo $P > ./tmp/disc.tmp


# background
awk '{s=0; d=0;  for(i=2;i<=NF;i++){s=s+$i; a[i]=$i} for(i=2;i<=NF;i++){d=d+(a[i]-s/(NF-1))^2;}  print sqrt(d/(NF-1))}' ./database/protein.dat | awk '{s=s+$1} END{print s/NR}' > ./tmp/varp.tmp
awk '{s=0; d=0;  for(i=2;i<=NF;i++){s=s+$i; a[i]=$i} for(i=2;i<=NF;i++){d=d+(a[i]-s/(NF-1))^2;}  print sqrt(d/(NF-1))}' ./database/rna.dat     | awk '{s=s+$1} END{print s/NR}' > ./tmp/varr.tmp
se2=`paste ./tmp/varp.tmp ./tmp/varr.tmp | awk 'BEGIN{s=0} ($1<0.005)||($2<0.02){s=1} END{print s}'`

if [ $se2 == 1 ]; then
V=0;
fi

# calculates the contributions
sh statistics.sh $1 $2
