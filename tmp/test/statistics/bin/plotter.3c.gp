#max=`sort -n -k 1 $1 | tail -1 | awk '{print int($1)+5}'`
#min=`sort -n -k 1 $1 | head -1 | awk '{print int($1)-5}'`
min=1; 
max=`awk '{print length($2)}' $5`
min1=`more $1 | sort -n -k 2 -r | awk '{print $2-0.5}' | tail -1`
max1=`more $1 | sort -n -k 2 | awk '{print $2+0.5}' | tail -1`



gnuplot << EOF
set style fill solid 0.50 
set border 3 lw 3 
set xtics nomirror
set ytics nomirror

set style fill pattern 3

# set terminal pdf 
# set output 'dp.pdf'

set terminal pngcairo font 'Arial, 14' size 900,562
set output 'dp.png'

set xlabel "$2 Position" font "Arial, 18"
set ylabel "Interaction Score"  font "Arial, 18"
set title "Interaction Profile" font "Arial Bold, 18"

set xrange [$min:$max] 
set yrange [$min1:$max1]


set xtics rotate by 0

#plot '$1' u 1:2 w histeps lc -1 notitle, 0 lc 1 lw 3 notitle
plot '$1' u 1:2  w filledcu y1=0 lc rgb "#B03B00"  notitle
EOF
