#rm *.pdf 
#!/bin/bash
PATH=$PATH:/sw/bin; export PATH
more ../download/interactions.$1.$2 | grep -v "#" | awk '($2!~/_R/){print $2}' | sort -k 1 | uniq > rna.names.txt
more ../download/interactions.$1.$2 | grep -v "#" | awk '($1!~/_R/){print $1}' | sort -k 1 | uniq > protein.names.txt

for j1 in `cat protein.names.txt`; do for j2 in `cat rna.names.txt`; 
do awk '($1=="'$j1'")&&($2=="'$j2'"){s=$3} ($1~/'$j1'/)&&($2~/'$j2'/)&&(s!=$3){K++ ; if(s>$3){k++}} END{print "'$j1' + '$j2' Interaction Specificity:", int(k/K*100)"%","cases:", K, "score:", s}' ../download/interactions.$1.$2;
done; done > ../download/statistics.$1.$2
for j1 in `cat protein.names.txt`; do for j2 in `cat rna.names.txt`; 
do awk '($1=="'$j1'")&&($2=="'$j2'"){s=$3} ($1=="'$j1'")&&($2~/'$j2'/)&&(s!=$3){K++ ; if(s>$3){k++}}               END{print "'$j1' + '$j2' RNA         Specificity:", int(k/K*100)"%","cases:", K, "  score:",s}' ../download/interactions.$1.$2;
done; done >> ../download/statistics.$1.$2
for j1 in `cat protein.names.txt`; do for j2 in `cat rna.names.txt`; 
do awk '($1=="'$j1'")&&($2=="'$j2'"){s=$3} ($2=="'$j2'")&&($1~/'$j1'/)&&(s!=$3){K++ ; if(s>$3){k++}}               END{print "'$j1' + '$j2' Protein     Specificity:", int(k/K*100)"%","cases:", K, "  score:",s}' ../download/interactions.$1.$2; 
done; done >> ../download/statistics.$1.$2


# total
for j1 in `cat protein.names.txt`; do 
for j2 in `cat rna.names.txt`; do 
awk '($1=="'$j1'")&&($2=="'$j2'"){s=$3} ($2~/'$j2'/)&&($1~/'$j1'/)&&(s!=$3){print $3}' ../download/interactions.$1.$2  | awk 'BEGIN{
min=1000; max=-1000} (min>$1){min=$1} (max<$1){max=$1} {a[NR]=$1} END{
num=100; delta=(max-min)/num; b[1]=min; for(i=2;i<=num+1;i++){b[i]=b[i-1]+delta} for(i=1;i<=num;i++){c=0; 
for(j=1;j<=NR;j++){if((a[j]>=b[i])&&(a[j]<b[i+1])){c++}} print (b[i]+b[i+1])/2,c}}' > histo.total.txt; 
echo "Specificity for $j1 and $j2 (vs 10000 interactions)" | sed 's/_1//g' | awk '{print $1, $2, substr($3,1,15), $4, substr($5,1,15), $6, $7,$8}' > title.txt
sh ./statistics/bin/plotter.gp ../download/statistics.$1.$2 histo.total.txt 2 1 $j1 $j2
mv dp.pdf total.$j1.$j2.pdf
gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.pdf *.pdf
done ;done
# keeps memory
mv all.pdf all.memo
rm *.pdf

# rna
for j1 in `cat protein.names.txt`; do 
for j2 in `cat rna.names.txt`; do 
awk '($1=="'$j1'")&&($2=="'$j2'"){s=$3} ($1=="'$j1'")&&($2~/'$j2'/)&&(s!=$3){print $3}' ../download/interactions.$1.$2 | awk 'BEGIN{
min=1000; max=-1000} (min>$1){min=$1} (max<$1){max=$1} {a[NR]=$1} END{
num=10; delta=(max-min)/num; b[1]=min; for(i=2;i<=num+1;i++){b[i]=b[i-1]+delta} for(i=1;i<=num;i++){c=0; 
for(j=1;j<=NR;j++){if((a[j]>=b[i])&&(a[j]<b[i+1])){c++}} print (b[i]+b[i+1])/2,c}}' > histo.rna.txt 
echo "$j2 Specificity for $j1 (vs 100 RNAs)" | sed 's/_1//g' | awk '{print substr($1,1,15), $2, $3, substr($4,1,15), $5,$6, $7}' > title.txt
sh ./statistics/bin/plotter.gp ../download/statistics.$1.$2 histo.rna.txt 50 2 $j1 $j2
mv dp.pdf rna.$j1.$j2.pdf
gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.pdf *.pdf 
done ; done 
# keeps memory
rm rna*.pdf
mv all.pdf  all.2.pdf
mv all.memo all.1.pdf
gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.pdf *.pdf 
mv all.pdf all.memo
rm *.pdf

# protein
for j1 in `cat protein.names.txt`; do 
for j2 in `cat rna.names.txt`; do 
awk '($1=="'$j1'")&&($2=="'$j2'"){s=$3} ($2=="'$j2'")&&($1~/'$j1'/)&&(s!=$3){print $3}' ../download/interactions.$1.$2 | awk 'BEGIN{
min=1000; max=-1000} (min>$1){min=$1} (max<$1){max=$1} {a[NR]=$1} END{
num=10; delta=(max-min)/num; b[1]=min; for(i=2;i<=num+1;i++){b[i]=b[i-1]+delta} for(i=1;i<=num;i++){c=0; 
for(j=1;j<=NR;j++){if((a[j]>=b[i])&&(a[j]<b[i+1])){c++}} print (b[i]+b[i+1])/2,c}}' > histo.protein.txt
echo "$j1 Specificity for $j2 (vs 100 proteins)"| sed 's/_1//g'| awk '{print substr($1,1,15), $2, $3, substr($4,1,15), $5,$6, $7}' > title.txt
sh ./statistics/bin/plotter.gp ../download/statistics.$1.$2 histo.protein.txt 50 3 $j1 $j2
mv dp.pdf protein.$j1.$j2.pdf
gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.pdf *.pdf
done ; done 
rm protein*.pdf
## keeps memory
mv all.pdf  all.2.pdf
mv all.memo all.1.pdf
gs -q -sPAPERSIZE=a4 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.pdf *.pdf 
mv all.pdf ../download/report.$1.$2.pdf
rm  *.pdf
