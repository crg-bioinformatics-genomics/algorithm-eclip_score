#!/bin/bash
mkdir ./tmp
mkdir ./results
mkdir ./database

name=`echo $1 | sed 's/\.txt//g' | awk '{print $1}'`

if [[ -s "$name.rna.frag.lib" ]]; then
	rm "$name.rna.frag.lib"
fi

if [[ -s "$name.rna.lib" ]]; then
	rm "$name.rna.lib"
fi

if [[ -s "errors.log" ]]; then
	rm "errors.log"
fi

# run cases one by one
for i2 in `cat $1 | grep -v "#" | awk '{print $1}' | head -20000 | sed 's/>//g'`; do 

	si2=`grep -w $i2 $1 | awk '(NF==2){print $2}'`
		
	if [[ ${#si2} -lt 50 ]]; then 

		echo "$i2" >> errors.log

	elif [[ ${#si2} -gt 1200 ]]; then 
	
		bash ./bin/start.fragment.sh $i2 $si2 >> $name.rna.frag.lib
		
		cat tmp/rna.frags >> $name.frag.txt
	
	else 
		
		bash ./bin/start.modified.sh $i2 $si2 >> $name.rna.lib

	fi
	
done

rm -fr tmp
rm -fr database
rm -fr results
