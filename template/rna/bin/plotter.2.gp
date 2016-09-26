max=`wc $1 | awk '{print $1+1}'`
fon=`wc $1 | awk '{print 2.5+(1-$1/100)*2.5}'`
gnuplot << EOF
set style fill solid 0.50 noborder

set terminal pdf 
# nocrop enhanced  font arial size 1024,768
set output 'dp.pdf'

set xlabel "$2 Position" 
set ylabel "Interaction Score"
#set title '$2'
#set parametric
#const=$start
#set trange [0:$height]
set xrange [0:$max] 
set xtics rotate by -45 font "Helvetica, $fon"
#set yrange [0:$height]
#set label " Score = $score  -   Significance = $perc"        at first  $start+1, first $height/3+0.2 front rotate by 90
#set arrow from $start,$height/3+0.2 to $length+$diff,$height/3+0.2   head front nofilled lc -1  linewidth 3
plot '$1' u 1:3:4:xtic(2) w errorline lc -1 lw $3 notitle 
EOF
