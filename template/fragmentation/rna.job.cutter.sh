file=$1

for i in `awk '{print $1}' $file | sort -u`; do
	awk '($1=="'$i'")' $file > tmp.txt
	nt=`cat tmp.txt  | awk '(NR==1){s=length($2)/50; if(s<=25){s=25}   if(s>=750){s=750} print int(s)}'`
	bash cutter.sh tmp.txt $nt | sort -u | awk '{print "'$i'_"$1,$2}'
	rm tmp/cuts.tmp.txt.*
done
