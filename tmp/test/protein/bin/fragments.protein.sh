#s=$2; 
#echo "delete" > tmp/outy.txt
#echo "delete" > tmp/outz.txt
#c=0;
 
#while [ $c == 0 ] && [ $s != 50 ]; 

#do
#more $1 | awk '{for(i=1;i<=length($2);i=i+int('$s'/7)){ print i, substr($2,i,'$s')}}' | awk '(length($2)>='$s'){ch=0; for(i=1;i<=length($2);i++){aa=substr($2,i,1);  if(aa=="K"){ch++}
#if(aa=="R"){ch++} if(aa=="E"){ch--} if(aa=="D"){ch--}} end='$s'+$1; print $1"-"end, ch/length($2), $2}' | sort -n -k 2 | uniq | awk '($2>=0){print $1, $3}' > tmp/outy.txt

#b=`wc tmp/outy.txt | awk '($1>100){print 1} ($1<=100){print 0}'`
#if [ $b == 0 ] ; then 
#cat tmp/outy.txt > tmp/outz.txt
#fi 

#c=`wc tmp/outz.txt | awk '($1>100){print 1} ($1<=100){print 0}'`
#let s=s-10;

#done
#grep -v "delete" tmp/outz.txt



s=50; 
echo "delete" > tmp/outy.txt
echo "delete" > tmp/outz.txt
c=0; 
while [ $c == 0 ] && [ $s != 200 ]; 

do
more $1 | awk '{for(i=1;i<=length($2);i=i+int('$s'/5)){ print i, substr($2,i,'$s')}}' | awk '(length($2)>='$s'){ch=0; for(i=1;i<=length($2);i++){aa=substr($2,i,1);  if(aa=="K"){ch++}
if(aa=="R"){ch++} if(aa=="E"){ch--} if(aa=="D"){ch--}} end='$s'+$1; print $1"-"end, ch/length($2), $2}' | sort -n -k 2 | uniq | awk '($2>0){print $1, $3}' > tmp/outy.txt

b=`wc tmp/outz.txt | awk '($1>50){print 1} ($1<=50){print 0}'`
if [ $b == 0 ] ; then 
cat tmp/outy.txt >> tmp/outz.txt
fi 

c=`wc tmp/outz.txt | awk '($1>100){print 1} ($1<=100){print 0}'`
let s=s+10;

done
grep -v "delete" tmp/outz.txt

