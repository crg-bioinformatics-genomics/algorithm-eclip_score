
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

# dp 
sh ./bin/score.3.sh
awk '(NR==1){for(i=2;i<=NF;i++){a[i]=$i}} (NR>1){for(i=2;i<=NF;i++){s=s+$1*a[i]*$i}} END{print s}' ./tmp/rlc.tmp > ./tmp/score.2.txt
spec=`sh ./bin/processer.sh tmp/score.2.txt | awk '{print s=int($3*100)}'`

# red and blue contrast
params_path=./params
R=`more ./data/ff.$1-$2.mat.txt | awk '{for(i=1;i<=NF;i++){if($i>0){p=p+$i;P++} if($i<0){n=n+$i;N++}}} END{print p}'`
B=`more ./data/ff.$1-$2.mat.txt | awk '{for(i=1;i<=NF;i++){if($i>0){p=p+$i;P++} if($i<0){n=n+$i;N++}}} END{print n}'`

# ip
head -1     ./tmp/score.txt      > ./tmp/analysis.txt
V=`cat tmp/analysis.txt`
cat $params_path/positives.txt    >> ./tmp/analysis.txt
P=`awk '(NR==1){s=$1} (NR>1){if($1>s){k++}} END{print k/(NR-1)}' tmp/analysis.txt`

awk '('$V'>$1){k++} END{print k/NR}' $params_path/negatives.txt | awk '{

a=int(100*$1*(1-'$P')*(0.25)/((1-$1)*'$P'*(0.75)+(1-'$P')*$1*(0.25)))
b=int(100*$1*(1-'$P')*(0.50)/((1-$1)*'$P'*(0.50)+(1-'$P')*$1*(0.50)))
c=int(100*$1*(1-'$P')*(0.75)/((1-$1)*'$P'*(0.25)+(1-'$P')*$1*(0.75)))
d=int(100*$1*(1-'$P')*(0.88)/((1-$1)*'$P'*(0.12)+(1-'$P')*$1*(0.88)))
e=int(100*$1*(1-'$P')*(0.33)/((1-$1)*'$P'*(0.67)+(1-'$P')*$1*(0.33)))

print "'$1'+'$2'", "Propensity:", a"%", e"%", b"%", c"%", d"%", "Red to Blue Ratio:", '$R', '$B'}' | grep -v "#" | sed 's/%//g' | awk '{
#
slope=0.06;
intercept=-3.44;
fslope=98.2;
fintercept=48.3;
v=$4;
#
r=(slope*$(NF+0)*(+1)+intercept/2);  if(r<0){sr=-int(r*100)/100 } if(r>=0){sr=int(r*100)/100}  
b=(slope*$(NF-1)*(+1)+intercept/2);  if(b<0){sb=-int(b*100)/100 } if(b>=0){sb=int(b*100)/100}
#
if(sr>sb){arb=sr} if(sr<=sb){arb=sb}
r2=0; b2=0;
#
if(r<0){r2=r+arb} if(r>=0){r2=r-arb}
if(b<0){b2=b+arb} if(b>=0){b2=b-arb}
#
r=r2;
b=b2;
rb=(slope*($NF+$(NF-1))+intercept);  tbr=(exp(rb)-exp(-rb))/(exp(rb)+exp(-rb));  
#
tr=(exp(r)-exp(-r))/(exp(r)+exp(-r)); tb=(exp(b)-exp(-b))/(exp(b)+exp(-b)); 
B=(tr)/(1+tr*tb); R=(tb)/(1+tr*tb);
#
SR=(R+1)/2; 
SB=(B-1)/2;
#
SR=(SR)*fslope+fintercept/2;
SB=(SB)*fslope+fintercept/2;
#
if(SR>0){d1=0}  if(SR<=0){d1=SR} 
if(SB>0){d2=SB} if(SB<=0){d2=0} 
#
SR=SR-d1-d2
SB=SB+d1-d2
#
#if(SB==0){SR=v}
#t=(SR+SB)/v;

rem=v-(SR+SB)
if(rem>0){SR=SR+rem}
if(rem<0){SB=SB+rem}

#SR=+SR/t;
#SB=v-SR
 
if(SB<0.001){SB=-SB}

#changed
SR=int('$spec');
if((v<=50)&&(SR<=50)){sc1="<  33%";   sc2="<  33%";}
if((v>50 )&&(SR> 50)){sc1=" = "v"%";  sc2=" = "SR"%";}
if((v<=50)&&(SR>=50)){sc1="<= 50% -"v"%";   sc2=" = "SR"%";}
if((v>50 )&&(SR< 50)){sc1=" = "v"%";  sc2="<= 50% -"SR"%";}

print "'$1'+'$2'", "Interaction Propensity ", sc1, "Discriminative Power ", sc2}' 
