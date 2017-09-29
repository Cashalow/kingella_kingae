#!/bin/bash

#first argument is prot or gene, second argument is raxml model

type=$1
model=$2
target_folder=$3

index=0
start=1
rm raxml_partitionfile.txt
for i in $(ls ${target_folder}/*_oneline)
do
  length=$(tail -n -1 ${i} | tr -d '\n' | wc -m);
  echo "${model}, ${type}${index} = ${start}-$(( start + length - 1 ))" >> ${target_folder}/raxml_partitionfile.txt
  start=$(( start + length ));
  index=$(( index + 1 ));
done


rm ${target_folder}/*_oneline
