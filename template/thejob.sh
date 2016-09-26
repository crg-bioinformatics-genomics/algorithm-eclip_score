for i in `ls ../rna_seqs/`;
do
  sed 's/[\| | \\ | \/ | \* | \$ | \# | \? | \! ]/./g' ../rna_seqs/$i | awk '(length($1)>=1)' | awk '($1~/>/){gsub(" ", "."); printf "\n%s\t", $1} ($1!~/>/){gsub(/[Uu]/, "T", $1); printf "%s",toupper($1)}' | awk '(NF==2)' | head -1 | sed 's/>\.//g;s/>//g' | awk '{print substr($1,1,12)"_"NR, $2}' >  rna_seqs_oneline/$i
done

cd fragmentation
for i in `ls ../rna_seqs_oneline/`;
do
  bash rna.job.cutter.sh ../rna_seqs_oneline/$i > ../rna_seqs_fragmented/$i.rna.fragm.seq;
done
cd ..

cd rna_seqs_fragmented
n=0
for i in *
do
  if [ $((n+=1)) -gt 10 ]; then
    n=1
  fi
  todir=../rna.libraries.U/rna_seqs/$n
  [ -d "$todir" ] || mkdir "$todir" 
  mv "$i" "$todir" 
done

cd ../rna.libraries.U/
bash superjob.sh
cd ../
