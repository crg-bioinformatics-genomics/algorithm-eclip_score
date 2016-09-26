score=`more $1 | grep -w $5 | grep -w $6 | tail -1  | awk '{print $NF}'`
perc=`more $1 | grep -w $5 | grep -w $6 | head -$4 | tail -1 | awk '{print $6}'`
height=`more $2 | sort -n -k 2 | tail -1 | awk '{print $NF}'`
minhe=`more $2 | sort -n -k 1 | head -1 | awk '{print $1}'`
length=`more $2 | sort -n -k 1 | tail -1 | awk '($1>'$score'){print $1} ($1<'$score'){print '$score'}'`
diff=`echo $length $minhe | awk '{print ($1-$2)/2}'`
title=`cat title.txt`  
start=`echo $score $minhe | awk '($1<$2){print $2} ($1>=$2){print $1}'`

gnuplot << EOF
set style fill solid 0.50 noborder

set terminal pdf 
# nocrop enhanced  font arial size 1024,768
set output 'dp.pdf'

set xlabel "Interaction Score"
set ylabel "Number of Counts"
set title '$title : Score = $score  /  Significance = $perc'
set parametric
const=$start
set trange [0:$height]
set xrange [$minhe-$diff:$length+$diff]
set yrange [0:$height]
#set label " Score = $score  -   Significance = $perc"        at first  $start+1, first $height/3+0.2 front 

set arrow from $start,$height/3+0.2 to $length+$diff,$height/3+0.2   head front nofilled lc 3  linewidth 3
plot const,t lc 3 with filledcurve x1=$length+$diff notitle, '$2' w i lc -1 lw $3 notitle
EOF
