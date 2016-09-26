echo "# Protein Sequences" > ../tmp/prot.names.txt
echo "# RNA Sequences" > ../tmp/rna.names.txt
echo " " > ../tmp/space.tmp
cat ../tmp/prot.names.txt ../tmp/space.tmp ../protein/outfile.fr ../tmp/space.tmp  ../tmp/rna.names.txt ../tmp/space.tmp  ../rna/outfile2.fr >> ../outputs/yoursequences.$1.$2.txt

grep -v "#" ../outputs/interactions.$1.$2.txt | grep -v "#"| awk '{print $3, $2}' | sed 's/_/ /g' | awk '{print $NF, $1}'  | sort -n -k 1 > ../tmp/rna.tmp
grep -v "#" ../outputs/interactions.$1.$2.txt | grep -v "#"| awk '{print $3, $1}' | sed 's/_/ /g' | awk '{print $NF, $1}'  | sort -n -k 1 > ../tmp/prot.tmp

bash bin/norm.2.sh ../tmp/prot.tmp > ../tmp/prot.2.tmp
bash bin/sorter.sh ../tmp/prot.tmp > ../tmp/prot.s.tmp

cat ../outputs/interactions.$1.$2.txt | awk '(NR==1){print $0} (NR>1){s=s+$3;a[NR]=$3; b[NR]=$1;  c[NR]=$2; e[NR]=$4}
END{for(i=2;i<=NR;i++){d=d+(a[i]-s/NR)^2} for(i=2;i<=NR;i++){print b[i], c[i], (a[i]-s/NR)/sqrt(d/NR), e[i],a[i]}}' | grep -v "#" | sort -n -k 3 | tail -20 | awk '{print NR, $1"+"$2, $3,$5}' > ../tmp/best.prot.rna.tmp
cat ../outputs/interactions.$1.$2.txt | awk '(NR==1){print $0} (NR>1){s=s+$3;a[NR]=$3; b[NR]=$1;  c[NR]=$2; e[NR]=$4}
END{for(i=2;i<=NR;i++){d=d+(a[i]-s/NR)^2} for(i=2;i<=NR;i++){print b[i], c[i], (a[i]-s/NR)/sqrt(d/NR), e[i],a[i]}}' | grep -v "#" | sort -n -k 3 | tail -20 | sort -nrk3 | awk '{print NR, $1,$2,$5,$4,$3}' > ../outputs/best.prot.rna.txt

gs=`awk '{print ($NF+1)/2}' ../filter/processed.txt`
awk '{print $1,$2,$3,$4,$5,$6+'$gs'}' ../outputs/best.prot.rna.txt > ../outputs/best.prot.rna.txt2
mv ../outputs/best.prot.rna.txt2 ../outputs/best.prot.rna.txt 

cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2.txt
echo "# Interactions - Raw Scores " >> ../outputs/yoursequences.$1.$2.txt
cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2.txt
cat  ../tmp/best.prot.rna.tmp | awk '{print $2, $4}' >> ../outputs/yoursequences.$1.$2.txt


cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2.txt
echo "# Protein - Raw Scores " >> ../outputs/yoursequences.$1.$2.txt
cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2.txt
cat  ../tmp/prot.s.tmp >> ../outputs/yoursequences.$1.$2.txt

bash bin/norm.2.sh ../tmp/rna.tmp >  ../tmp/rna.2.tmp
bash bin/sorter.sh ../tmp/rna.tmp >  ../tmp/rna.s.tmp

cat ../tmp/rna.2.tmp | awk '{print $1,$2*'$gs'}' > ../tmp/rna.2.tmp2
mv ../tmp/rna.2.tmp2 ../tmp/rna.2.tmp
cat ../tmp/rna.s.tmp | awk '{print $1,$2*'$gs'}' > ../tmp/rna.s.tmp2
mv ../tmp/rna.s.tmp2 ../tmp/rna.s.tmp

cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2.txt
echo "# RNA - Raw Scores" >> ../outputs/yoursequences.$1.$2.txt
cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2.txt
cat  ../tmp/rna.s.tmp >> ../outputs/yoursequences.$1.$2.txt

cp ../outputs/interactions.$1.$2.txt int.tmp
awk '(NR==1){print $0} (NR>1){print $1,$2,$3*'$gs',$4}' int.tmp > int.tmp2
mv int.tmp2 int.tmp

bash bin/map.sh int.tmp $1 $2

# mv dp.pdf 5best.pdf

# plots 
# bash bin/plotter.3.gp  ../tmp/best.prot.rna.tmp Score 3 Best
# mv dp.pdf 1best.pdf
bash bin/plotter.3c.gp ../tmp/rna.2.tmp   Nucleotide 3 RNA ../rna/outfile2
mv dp.png ../outputs/Interaction_profile.$1.$2.png
# mv dp.pdf 3rna.pdf
# bash bin/plotter.3c.gp ../tmp/prot.2.tmp  Aminoacid  3 Protein ../protein/outfile
# mv dp.pdf 4protein.pdf
# gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.pdf *.pdf 
# cp all.pdf ../outputs/report.$1.$2.pdf
