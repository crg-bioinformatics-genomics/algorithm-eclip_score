n=10
for i in `ls rna_seqs/$n/*.seq | xargs -n1 basename`; do
cd rna_library_generator$n/
sed 's/|/_/g' ../rna_seqs/$n/"$i" > "$i".txt
bash rungenerator.rna.sh "$i".txt
mv "$i".rna.lib ../outs
rm "$i".txt
cd ..
done
