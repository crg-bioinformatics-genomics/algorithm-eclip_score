echo "# Protein Sequences" > ../tmp/prot.names.txt
# PATH=$PATH:/sw/bin; export PATH
echo "# RNA Sequences" > ../tmp/rna.names.txt
echo " " > ../tmp/space.tmp
cat ../tmp/prot.names.txt ../tmp/space.tmp ../protein/outfile.fr ../tmp/space.tmp  ../tmp/rna.names.txt ../tmp/space.tmp  ../rna/outfile2.fr > ../outputs/yoursequences.$1.$2

grep -v "#" ../outputs/interactions.$1.$2.txt | grep -v "#"| awk '{print $3, $2}' | sed 's/_/ /g' | awk '{print $NF, $1}'  | sort -n -k 1 > ../tmp/rna.tmp
grep -v "#" ../outputs/interactions.$1.$2.txt | grep -v "#"| awk '{print $3, $1}' | sed 's/_/ /g' | awk '{print $NF, $1}'  | sort -n -k 1 > ../tmp/prot.tmp

for j in `cat ../tmp/prot.tmp | awk '{print $1}' | uniq`; do 
#awk '($1=="'$j'")' ../tmp/prot.tmp | awk '{a[NR]=$2; s=s+$2} END{for(i=1;i<=NR;i++){d=d+(a[i]-s/NR)^2} print "'$j'", s/NR, sqrt(d)/NR}'; 
awk '($1=="'$j'")' ../tmp/prot.tmp | sort -n -k 2 | tail -1 | awk '{ print "'$j'", $2}'
done | sort -n -k 2 | awk '{print NR, $0}' > ../tmp/prot.2.tmp

grep -v "#" ../outputs/interactions.$1.$2.txt | sort -n -k 3    | tail -20 | awk '{print NR, $1"+"$2, $3}' > ../tmp/best.prot.rna.tmp
cat ../outputs/interactions.$1.$2.txt | awk '(NR==1){print $0} (NR>1){s=s+$3;a[NR]=$3; b[NR]=$1;  c[NR]=$2; e[NR]=$4}END{for(i=2;i<=NR;i++){d=d+(a[i]-s/NR)^2} for(i=2;i<=NR;i++){print b[i], c[i], (a[i]-s/NR)/sqrt(d/NR), e[i],a[i]}}' | grep -v "#" | sort -n -k 3 | tail -20 | sort -nrk3 | awk '{print NR, $1,$2,$5,$4,$3}' > ../outputs/best.prot.rna.txt
grep -v "#" ../outputs/interactions.$1.$2.txt | sort -n -k 3 -r | tail -20 | awk '{print NR, $1"+"$2, $3}' > ../tmp/worst.prot.rna.tmp

gs=`awk '{print ($NF+1)/2}' ../filter/processed.txt`
awk '{print $1,$2,$3,$4,$5,$6+'$gs'}' ../outputs/best.prot.rna.txt > ../outputs/best.prot.rna.txt2
mv ../outputs/best.prot.rna.txt2 ../outputs/best.prot.rna.txt 

cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2
echo "# Interactions - Best Scores " >> ../outputs/yoursequences.$1.$2
cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2
cat  ../tmp/best.prot.rna.tmp >> ../outputs/yoursequences.$1.$2


cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2
echo "# Protein - Max Scores " >> ../outputs/yoursequences.$1.$2
cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2
cat  ../tmp/prot.2.tmp >> ../outputs/yoursequences.$1.$2

for j in `cat ../tmp/rna.tmp  | awk '{print $1}' | uniq`; do 
#awk '($1=="'$j'")' ../tmp/rna.tmp | awk '{a[NR]=$2; s=s+$2} END{for(i=1;i<=NR;i++){d=d+(a[i]-s/NR)^2} print "'$j'", s/NR, sqrt(d)/NR}'; 
awk '($1=="'$j'")' ../tmp/rna.tmp | sort -n -k 2 | tail -1 | awk '{ print "'$j'", $2}'
done  | sort -n -k 2 | awk '{print NR, $0}' > ../tmp/rna.2.tmp

cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2
echo "# RNA - Max Scores" >> ../outputs/yoursequences.$1.$2
cat ../tmp/space.tmp >> ../outputs/yoursequences.$1.$2
cat  ../tmp/rna.2.tmp >> ../outputs/yoursequences.$1.$2

cat  ../outputs/interactions.$1.$2.txt | sed 's/-/ -/g'  | grep -v "#" | awk '(-$2-$3>120){print $(NF-3), $(NF-1)}' > ../tmp/data.txt

cat ../tmp/data.txt | sort -n -k 1 | awk 'BEGIN{print 0, 0, 0} {a[NR]=$1; z[NR]=$2} 
END{for(i=1;i<=NR;i++){

#if(i==1){print a[i],  z[i]}
#if(i==2){print a[i], (z[i]+z[i-1]+z[i+1])/3}
#if(i==3){print a[i], (z[i]+z[i-1]+z[i+1]+z[i+2]+z[i-2])/5}
#if(i==NR)  {print a[i], z[i]}
#if(i==NR-1){print a[i], (z[i]+z[i-1]+z[i+1])/3}
#if(i==NR-2){print a[i], (z[i]+z[i-1]+z[i+1]+z[i+2]+z[i-2])/5}

if((i>3)&&(i<NR-2)){print (a[i-3]+a[i-2]+a[i-1]+a[i]+a[i+1]+a[i+2]+a[i+3])/7, (z[i-3]+z[i-2]+z[i-1]+z[i]+z[i+1]+z[i+2]+z[i+3])/7,0}
}}' ../tmp/data.txt  | sort -n -k 1 > ../tmp/smoothed.txt

awk '{s=s+$2; a[NR]=$1; b[NR]=$2} END{for(i=1;i<=NR;i++){d=d+(b[i]-s/NR)^2} for(i=1;i<=NR;i++){print a[i], (b[i]-s/NR)/sqrt(d/NR),b[i]}}'  ../tmp/smoothed.txt  | awk '($3>=5.03)||(($3<5.03)&&($2<0)){print $0,0} ($3<5.03)&&($2>0){print $1, 0, $3,0}'> ../tmp/smoothed.2.txt

cat ../tmp/smoothed.txt | awk '{print $1,$2*'$gs',$3,$4}' > ../tmp/smoothed.txt2
mv ../tmp/smoothed.txt2 ../tmp/smoothed.txt
cat ../tmp/smoothed.2.txt | awk '{print $1,$2*'$gs',$3,$4}' > ../tmp/smoothed.2.txt2
mv ../tmp/smoothed.2.txt2 ../tmp/smoothed.2.txt

cp ../tmp/smoothed.2.txt ../outputs/smoothed.z.$1.$2
cp ../tmp/smoothed.txt ../outputs/smoothed.$1.$2

bash bin/plotter.5.gp ../tmp/smoothed.txt
mv dp.png ../outputs/Interaction_profile.$1.$2.png
bash bin/plotter.6.gp ../tmp/smoothed.2.txt Z-
mv dp.png ../outputs/Interaction_Z-profile.$1.$2.png

# gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.pdf *.pdf 
# cp all.pdf ../outputs/report.$1.$2.pdf

#sh map.sh ../outputs/interactions.$1.$2.txt
#mv dp.pdf 5best.pdf
#
## plots 
#sh plotter.3.gp  ../tmp/best.prot.rna.tmp Score 3 Best
#mv dp.pdf 1best.pdf
#sh plotter.3.gp  ../tmp/worst.prot.rna.tmp Score 3 Worst
#mv dp.pdf 2worst.pdf
#sh plotter.3.gp ./tmp/rna.2.tmp   Nucleotide 3 RNA
#mv dp.pdf 3rna.pdf
#sh plotter.3.gp ./tmp/prot.2.tmp  Aminoacid  3 Protein
#mv dp.pdf 4protein.pdf
#gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.pdf *.pdf 
#cp all.pdf ../outputs/report.$1.$2.pdf
