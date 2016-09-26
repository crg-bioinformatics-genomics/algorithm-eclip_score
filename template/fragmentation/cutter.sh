awk '{for(i=1;i<=length($2);i=i+'$2'){if(i+'$2'<length($2)-'$2'){print $1, i,substr($2,i,'$2'*2+1)}}}'   $1 >  tmp/cuts.$1.$2
awk '{for(i=length($2);i>=1;i=i-'$2'){if(i-'$2'>='$2'+1){print $1, i-'$2'*2, substr($2,i-'$2'*2,'$2'*2+1)}}}' $1 >> tmp/cuts.$1.$2
cat tmp/cuts.$1.$2 |  sort -n -k 2 | awk '{print $2"-"'$2'*2+1+$2, $3}' 
