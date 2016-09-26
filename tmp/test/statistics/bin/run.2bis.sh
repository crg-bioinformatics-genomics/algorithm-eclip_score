echo "# Protein Sequences" > tmp/prot.names.txt
PATH=$PATH:/sw/bin; export PATH
echo "# RNA Sequences" > tmp/rna.names.txt
echo " " > tmp/space.tmp
cat tmp/prot.names.txt tmp/space.tmp outfile.fr tmp/space.tmp  tmp/rna.names.txt tmp/space.tmp  outfile2.fr >> ../download/yoursequences.$1.$2

grep -v "#" ../download/interactions.$1.$2 | grep -v "#"| awk '{print $3, $2}' | sed 's/_/ /g' | awk '{print $NF, $1}'  | sort -n -k 1 > tmp/rna.tmp
grep -v "#" ../download/interactions.$1.$2 | grep -v "#"| awk '{print $3, $1}' | sed 's/_/ /g' | awk '{print $NF, $1}'  | sort -n -k 1 > tmp/prot.tmp

for j in `cat tmp/prot.tmp | awk '{print $1}' | uniq`; do 
#awk '($1=="'$j'")' tmp/prot.tmp | awk '{a[NR]=$2; s=s+$2} END{for(i=1;i<=NR;i++){d=d+(a[i]-s/NR)^2} print "'$j'", s/NR, sqrt(d)/NR}'; 
awk '($1=="'$j'")' tmp/prot.tmp | sort -n -k 2 | tail -1 | awk '{ print "'$j'", $2}'
done | sort -n -k 2 |  awk '{a[NR]=$0} END{for(i=1;i<=25;i++){print a[i]} for(i=NR-24;i<=NR;i++){print a[i]}}' |  sort -n -k 2 | uniq |  awk '{print NR, $0}'  > tmp/prot.2.tmp

grep -v "#" ../download/interactions.$1.$2 | sort -n -k 3    | tail -20 | awk '{print NR, $1"+"$2, $3}' > tmp/best.prot.rna.tmp
grep -v "#" ../download/interactions.$1.$2 | sort -n -k 3 -r | tail -20 | awk '{print NR, $1"+"$2, $3}' > tmp/worst.prot.rna.tmp

cat tmp/space.tmp >> ../download/yoursequences.$1.$2
echo "# Interactions - Best Scores " >> ../download/yoursequences.$1.$2
cat tmp/space.tmp >> ../download/yoursequences.$1.$2
cat  tmp/best.prot.rna.tmp >> ../download/yoursequences.$1.$2


cat tmp/space.tmp >> ../download/yoursequences.$1.$2
echo "# Protein - Max Scores " >> ../download/yoursequences.$1.$2
cat tmp/space.tmp >> ../download/yoursequences.$1.$2
cat  tmp/prot.2.tmp >> ../download/yoursequences.$1.$2

for j in `cat tmp/rna.tmp  | awk '{print $1}' | uniq`; do 
#awk '($1=="'$j'")' tmp/rna.tmp | awk '{a[NR]=$2; s=s+$2} END{for(i=1;i<=NR;i++){d=d+(a[i]-s/NR)^2} print "'$j'", s/NR, sqrt(d)/NR}'; 
awk '($1=="'$j'")' tmp/rna.tmp | sort -n -k 2 | tail -1 | awk '{ print "'$j'", $2}'
done  | sort -n -k 2 | awk '{a[NR]=$0} END{for(i=1;i<=25;i++){print a[i]} for(i=NR-24;i<=NR;i++){print a[i]}}' | sort -n -k 2 | uniq | awk '{print NR, $0}' > tmp/rna.2.tmp

cat tmp/space.tmp >> ../download/yoursequences.$1.$2
echo "# RNA - Max Scores" >> ../download/yoursequences.$1.$2
cat tmp/space.tmp >> ../download/yoursequences.$1.$2
cat  tmp/rna.2.tmp >> ../download/yoursequences.$1.$2

sh map.sh ../download/interactions.$1.$2
mv dp.pdf 5best.pdf

# plots 
sh plotter.3.gp  tmp/best.prot.rna.tmp Score 3 Best
mv dp.pdf 1best.pdf
sh plotter.3.gp  tmp/worst.prot.rna.tmp Score 3 Worst
mv dp.pdf 2worst.pdf
sh plotter.3b.gp ./tmp/rna.2.tmp   Nucleotide 3 RNA
mv dp.pdf 3rna.pdf
sh plotter.3b.gp ./tmp/prot.2.tmp  Aminoacid  3 Protein
mv dp.pdf 4protein.pdf
gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.pdf *.pdf 
cp all.pdf ../download/report.$1.$2.pdf
