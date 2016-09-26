file=$1

for i in `awk '{print $1}' $file | sort -u`; do
	awk '($1=="'$i'")' $file > tmp.txt
	pr=`cat tmp.txt  | awk '(NR==1){s=length($2)/50; if(s<=25){s=25} if(s>=375){s=375} print int(s)}'`
	bash cutter.sh tmp.txt $pr | sort -u | awk '{print "'$i'_"$1,$2}'
	rm tmp/cuts.tmp.txt.*
done