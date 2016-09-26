min=`more ../tmp/smoothed.2.txt | sort -n -k 2 -r | awk '{print $2-0.5}' | tail -1`
max=`more ../tmp/smoothed.2.txt | sort -n -k 2 | awk '{print $2+0.5}' | tail -1`


gnuplot << EOF
set style fill solid 0.50 
#set style data lines
set style fill pattern 3
# set terminal pdf 
# set output 'dp.pdf'

set terminal pngcairo font 'Arial, 14' size 900,562
set output 'dp.png'


set border 3 lw 3 
set xtics nomirror
set ytics nomirror


set xlabel "Nucleotide Position" font "Arial, 18"
set ylabel "Interaction $2Score" font "Arial, 18"
set title "Interaction Z-Profile" font "Arial Bold, 18"
set yrange [$min:$max]
plot '$1' u 1:2  w filledcu y1=0 lc rgb "#B03B00"  notitle
#plot '$1' u 1:2 w steps lc -1 notitle
#plot '$1' u 1:2 smooth bezier lw 2 notitle 
EOF
