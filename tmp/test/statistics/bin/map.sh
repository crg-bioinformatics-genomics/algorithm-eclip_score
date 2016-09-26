more $1  | awk '{print $1, $2}' > ../tmp/one.txt
sed 's/-/ /g' ../tmp/one.txt | awk '{print $1+($2-$1)/2, $3+($4-$3)/2}' > ../tmp/two.txt 
paste ../tmp/two.txt  $1  | grep -v "#" |  awk '{print $1, $2, $5, $NF}' | sort -n -k 2 | sort -n -k 1 > ../tmp/three.2.txt

for j in `more ../tmp/three.2.txt | awk '{print $1}' | sort -n -k 1 | uniq`; do 
	awk '('$j'==$1)' ../tmp/three.2.txt > ../tmp/med.txt; 
	for i in `sort -n -k 2 ../tmp/med.txt | awk '{print $2}' | uniq`; do 
		awk '($2=="'$i'")' ../tmp/med.txt | sort -n -k 3 | tail -1 | awk '($4>0.25){print $1, $2, $3} ($4<=0.25){print $1, $2, 0}'
	done | uniq
done > ../tmp/three.txt 
cat ../tmp/three.2.txt | sort -n -k 3 | tail -3 >> ../tmp/three.txt
cat ../tmp/three.txt | awk '{print $1"#"$2, $3}' | sort -k 1 | awk '{a[NR]=$1; b[NR]=$2} END{for(i=1;i<=NR;i++){if(a[i]!=a[i+1]){print a[i],b[i]}}}' | sed 's/#/ /g' | sort -n -k1 -k2 > ../tmp/three.tmp
cp ../tmp/three.tmp ../tmp/three.txt

sort -n -k 2 ../tmp/three.txt |awk '{a[NR]=$2} END{for(i=1;i<=NR;i++){if(a[i]!=a[i+1]){print a[i]}}}' | awk '{print NR, $1}' | awk '{a[NR]=$0} END{print a[1]; print a[NR]}' > ../tmp/x.txt
bash bin/correlator.bis.sh ../tmp/x.txt | grep "#" > ../tmp/coeff.x.txt
sort -n -k 1 ../tmp/three.txt |awk '{a[NR]=$1} END{for(i=1;i<=NR;i++){if(a[i]!=a[i+1]){print a[i]}}}' | awk '{print NR, $1}' | awk '{a[NR]=$0} END{print a[1]; print a[NR]}' > ../tmp/y.txt
bash bin/correlator.bis.sh ../tmp/y.txt | grep "#" > ../tmp/coeff.y.txt
awk '{a[NR]=$1; b[NR]=$3} END{for(i=1;i<=NR;i++){if(a[i]==a[i+1]){printf "%f\t",b[i]} if(a[i]!=a[i+1]){printf "%f\n",b[i]}}}' ../tmp/three.txt  > ../tmp/five.txt


gnuplot < ./bin/map.2.plt
mv dp.png ../outputs/Interaction_matrix.$2.$3.png