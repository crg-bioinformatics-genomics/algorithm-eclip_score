max=`wc $1 | awk '{print $1+2}'`
fon=`wc $1 | awk '{print 2+(1-$1/100)*2.5}'`
ymi=`sort -n -k 3 $1 | awk '{print $3-5}' | head -1`
yma=`sort -n -k 3 $1 | awk '{print $3+5}' | tail -1`

gnuplot << EOF
set style fill solid 0.50 
#noborder

set terminal pdf 
# nocrop enhanced  font arial size 1024,768
set output 'dp.pdf'

set xlabel "$2 Position" 
set ylabel "Interaction Score"
set title "$4 Scores"
set parametric
const=25.5;
set style line 1 lt 2 lc rgb "red" lw 3
set trange [$ymi:$yma]
set xrange [-1:$max] 
set xtics rotate by -45 font "Helvetica, $fon"
#set ytics rotate by -90 font "Helvetica, $fon"
#set yrange [$ymi:$yma]
#set label " Score = $score  -   Significance = $perc"        at first  $start+1, first $height/3+0.2 front rotate by 90
#set arrow from $start,$height/3+0.2 to $length+$diff,$height/3+0.2   head front nofilled lc -1  linewidth 3
plot const,t w d lc 6 lw 3 notitle, '$1' u 1:3:xtic(2) lt -1 notitle 
EOF
