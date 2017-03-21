

number_of_species=7
grep ">" *.fa -c | grep ".fa:${number_of_species}" | sed "s/:.*//" > aln_${number_of_species}.txt
rm aln_${number_of_species}_uniq.txt
for file_name in $(cat aln_${number_of_species}.txt)
do
  if [ `grep ">" ${file_name} | sed "s/_.*//" | uniq -d | wc -l` -eq 0 ]; then 
    echo "${file_name}" >> aln_${number_of_species}_uniq.txt
  else
     echo "${file_name}" has duplicates
  fi
  done;
mkdir raxml_files
cd raxml_files
for aln in $(cat ../aln_${number_of_species}_uniq.txt)
do
  echo ${aln}
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../${aln} > ${aln}_oneline;
  done;
sed -i "s/>\([0-9]\+\)_.*/>\1/" *
index=0
start=1
rm raxml_partitionfile.txt
for i in $(ls *_oneline)
do
  length=$(tail -n -1 ${i} | wc -m);
  echo "LG, prot${index} = ${start}-$(( start + length - 2 ))" >> raxml_partitionfile.txt
  start=$(( start + length - 1 ));
  index=$(( index + 1 ));
  done
  
paste -d'\0' *_oneline > core_genes_concatenated.fa
sed "s/>\([0-9]\+\)>.*/>\1/" sed "s/>\([0-9]\+\)>.*/>\1/"
