max=`wc $1 | awk '{print $1+2}'`
fon=`wc $1 | awk '{print 2+(1-$1/100)*2.5}'`
gnuplot << EOF
set style fill solid 0.50 

set terminal pdf 
set output 'dp.pdf'

set xlabel "$2 Position" 
set ylabel "Interaction Score"
set title "Top 20 Interactions"
set xrange [-1:$max] 
set xtics rotate by -45 font "Helvetica, $fon"
plot '$1' u 1:3:xtic(2) lt -1 notitle 
EOF
