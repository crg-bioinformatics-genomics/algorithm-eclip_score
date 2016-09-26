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
nr="u1861"
np="O75817"

#sr="GUUCGUGCUGAAGGCCUGUAUCCUAGGCUACACACUGAGGACUCUGUUCCUCCCCUUUCCGCCUAGGGGAAAGUCCCCGGACCUCGGGCAGAGAGUGCCACGUGCAUACGCACGUAGACAUUCCCCGCUUCCCACUCCAAAGUCCGCCAAGAAGCGUAUCCCGCUGAGCGGCGUGGCGCGGGGGCGUCAUCCGUCAGCUCCCUCUAGUUACGCAGGCAGUGCGUGUCCGCGCACCAACCACACGGGGCUCAUUCUCAGCGCGGCUG"
#sp="MAENREPRGAVEAELDPVEYTLRKRLPSRLPRRPNDIYVNMKTDFKAQLARCQKLLDGGARGQNACSEIYIHGLGLAINRAINIALQLQAGSFGSLQVAANTSTVELVDELEPETDTREPLTRIRNNSAIHIRVFRVTPK"
#
echo
rr=`echo $2 | tr a-z A-Z | sed 's/T/U/g' | awk '{for(i=1;i<=NF;i++){printf "%s",$i} } END{printf "\n"}'`
pp=`echo $1 | tr a-z A-Z | awk '{for(i=1;i<=NF;i++){printf "%s",$i} } END{printf "\n"}'`
echo "protein"
echo $pp | awk '{print substr($1,1,60)"..."}'
echo "rna"
echo $rr | awk '{print substr($1,1,60)"..."}'
echo 
echo
sh start.modified.sh $np $pp $nr $rr
