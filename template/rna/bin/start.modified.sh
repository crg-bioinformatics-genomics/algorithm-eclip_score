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
echo $2 | sed 's/T/U/g' | tr [a-z] [A-Z] | awk '{print "'$1'", $1}' > ../database/rna.50.3000.txt 
cp  ../database/rna.50.3000.txt  ../database/rna.txt
  
# computes 10 features (concatenated) secondary + polarity + hydro for rna and proteins
bash bin/rna-feature.sh        ../database/rna.50.3000.txt       > ../database/rna.dat

# normalizes the lengths
bash bin/adaptator.sh           ../database/rna.dat 3001         > ../database/rna.3000.dat

# computes fouriers' coefficients
bash bin/fourier.line.fast.sh    ../database/rna.3000.dat   | awk '(NF==53)&&($2!=0)'   

# > ./database/fourier.rna.3000.dat


