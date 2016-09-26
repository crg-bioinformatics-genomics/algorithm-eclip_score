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


# builds the matrix
sh score.2.sh

# computes the product
awk '(NR==1){for(i=2;i<=NF;i++){a[i]=$i}} (NR>1){for(i=2;i<=NF;i++){s=s+$1*a[i]*$i}} END{print s}' ./tmp/rlc.tmp > ./tmp/score.txt

# visualize the regions 
awk '(NR==1){for(i=2;i<=NF;i++){a[i]=$i}} (NR>1){for(i=2;i<=NF;i++){printf "%f\t",$1*a[i]*$i} printf "\n"}' ./tmp/rlc.tmp  > ./tmp/mat.tmp
cp ./tmp/mat.tmp ./data/ff.$1-$2.mat.txt
