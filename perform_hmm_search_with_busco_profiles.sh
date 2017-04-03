alignmentsid_folder=$1
group=$2 #this is the group from which the orthogroups are extracted
hmm_folder=$3 #this is the file where the hmmm profiles are

outfolder=${alignmentsid_folder}/hmmsearch_${group}/ #results files will be put here


mkdir ${outfolder}

for i in $(ls ${alignmentsid_folder}/*.fa_oneline)
do
    sed "s/-//g" ${i} > ${outfolder}/$(basename ${i})_nogap;
done #removing the gaps from the alignments to perform the hmm search


cat ${hmm_folder}/*.hmm  > ${outfolder}/hmm_profiles_${group}.hmm #concatenating all hmm profiles

hmmpress ${outfolder}/hmm_profiles_${group}.hmm

for j in ${outfolder}/*_nogap
do
    b=$(basename ${j})
    hmmscan -E 1e-10 --tblout ${outfolder}/${b}_resfile.txt ${outfolder}/hmm_profiles_${group}.hmm ${j} > /dev/null
done #performs the hmmmsearch with evalue = 1e-10




