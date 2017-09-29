#!/bin/bash


number_of_species=$1 #total number of species in the orthofinder run
alignmentsid_folder=$2 #folder with the alignments and numbers as sole ids for species

outfolder=${alignmentsid_folder}/coregenome/

mkdir ${outfolder}

grep ">" ${alignmentsid_folder}/*.fa* -c | grep ":${number_of_species}" | sed "s/:.*//" | xargs -I % basename %  > ${outfolder}/aln_${number_of_species}.txt #finds all orthogroups with entry number = number of species

rm ${outfolder}/aln_${number_of_species}_uniq.txt
for file_name in $(cat ${outfolder}/aln_${number_of_species}.txt)
do
  if [ `grep ">" ${alignmentsid_folder}/${file_name} | sed "s/_.*//" | uniq -d | wc -l` -eq 0 ]; then 
    echo "$( basename $file_name)" >> ${outfolder}/aln_${number_of_species}_uniq.txt
  else
     echo "${file_name}" has duplicates
  fi
  done; #finds orthogroups with duplicates (not one entry per species) and excludes them


rm ${outfolder}/aln_${number_of_species}.txt
