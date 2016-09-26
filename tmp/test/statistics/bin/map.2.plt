# set terminal pdf 
# set output 'dp.pdf'

set terminal pngcairo font 'Arial, 12' size 900,562
set output 'dp.png'
set pm3d map 

sy=`awk '{print $6}' ../tmp/coeff.y.txt`
iy=`awk '{print $5}' ../tmp/coeff.y.txt`
sx=`awk '{print $6}' ../tmp/coeff.x.txt`
ix=`awk '{print $5}' ../tmp/coeff.x.txt`
mx=`grep -v "#" int.tmp | sort -n -k 3    | tail -1 | awk '{print int($3)+7}'`
my=`grep -v "#" int.tmp | sort -n -k 3 -r | tail -1 | awk '{print int($3)-7}'`
rl=`awk '{print length($2)}' ../rna/outfile2`
pl=`awk '{print length($2)}' ../protein/outfile`



set palette model RGB defined (0 "white", 100 "white", 160 "white", 200 "red")
set xrange  [0:rl]
set yrange  [0:pl]
set cbrange  [my:mx]
set colorbox horiz user origin 500, 500

#set cbrange [my:mx]
#set palette model RGB defined (0 "white", 100 "white", 160 "white", 200 "red")
#set cbrange [2:*]
#unset colorbox
#set colorbox horiz user origin 500, 500

set ylabel "Protein Residue Index" font "Arial, 14"
set xlabel "RNA Nucleotide Index" font "Arial, 14"
set title "Interaction Matrix" font "Arial Bold, 16"
#plot "../tmp/five.txt" using (($1*sx+ix>0)?($1*sx+ix):0):(($2*sy+iy>0)?($2*sy+iy):0.1):3 matrix w image notitle
#plot "../tmp/five.txt" using  ($1*sx+ix):($2*sy+iy):((($3-my)/(mx-my)>160/200)?($3):0) matrix w image notitle
 plot "../tmp/five.txt" using  ($1*sx+ix):($2*sy+iy):(($3>my+160/200*(mx-my))?($3):0) matrix w image notitle
#plot "../tmp/five.txt" using  ($1*sx+ix):($2*sy+iy):3 matrix w image notitle

