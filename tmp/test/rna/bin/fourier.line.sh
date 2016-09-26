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
awk '{
printf "%s\t",$1;
j=1; for(i=2;i<=NF;i++){s[j]=$i;j++}
#
s[1]=0;
s[j]=s[1];  
#
for(k=0;k<j-1;k++){
c=0;
for(i=0;i<j-1;i++){
c=c+sqrt(2/(j-1))*s[i+1]*cos( (3.14/(j-1))*(k+1/2)*(i+1/2))};
 printf "%4.2f\t",c
}
printf "\n"}' $1
#; done 
