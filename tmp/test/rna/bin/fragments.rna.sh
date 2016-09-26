#!/bin/bash
awk '{print $2}' $1 > ../tmp/out
label=`awk '{print $1}' $1 `

s=2002 ; c=0;

#while [ $c == 0  ] ; 

#do
RNALfold  -L $s < ../tmp/out | awk '{print $0, length($1)}' | awk '(NF>=4){print $(NF-2), $(NF-1), $NF}' | awk '($3>50)' > ../tmp/out.2
# reports the sequences
for u in `cat ../tmp/out.2 | awk '{print $2"-"$3}'`; do 
	pos=`echo $u | sed 's/-/ /g' | awk '{print $1}'`;
	length=`echo $u | sed 's/-/ /g' | awk '{print $2}'`; 
	awk '{end='$pos'+'$length'; print '$pos'"-"end,substr($1,"'$pos'", "'$length'")}' ../tmp/out; 
done > ../tmp/outx.txt

#c=`wc tmp/outx.txt | awk '($1>'$s'){print 0} ($1<='$s'){print 1}'`
#let s=s-10; 
#wc tmp/outx.txt
#done

cat ../tmp/outx.txt | awk '(length($2)<'$3'+50)&&(length($2)>'$3'-50)'


