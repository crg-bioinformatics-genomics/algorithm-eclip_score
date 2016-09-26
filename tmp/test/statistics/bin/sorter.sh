for j in `cat $1  | awk '{print $1}' | uniq`; 
do 
awk '($1=="'$j'"){print "'$j'", $2}' $1 | sort -n -k 2 | tail -1 | sed 's/-/ -/g' | awk '{print "'$j'",$3}'
#; 
done | sort -n -k 2
