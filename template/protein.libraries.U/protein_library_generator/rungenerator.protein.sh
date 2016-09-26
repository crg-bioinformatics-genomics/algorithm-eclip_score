#!/bin/bash
mkdir ./tmp
mkdir ./results
mkdir ./database

name=`echo $1 | sed 's/\.txt//g' | awk '{print $1}'`

if [[ -s "$name.prot.lib" ]]; then
	rm "$name.prot.lib"
fi

if [[ -s "errors.log" ]]; then
	rm "errors.log"
fi

# run cases one by one
for i2 in `cat $1 | grep -v "#" | awk '{print $1}' | head -200000 | sed 's/>//g'`; do 

	si2=`grep -w $i2 $1 | awk '(NF==2){print $2}'`

	bash ./bin/start.modified.sh $i2 $si2 >> $name.prot.lib

done

rm -fr tmp
rm -fr database
rm -fr results
