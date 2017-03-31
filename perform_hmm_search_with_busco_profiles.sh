alignmentsid_folder=$1
group=$2 #this is the group from which the orthogroups are extracted
hmm_folder=$3 #this is the file where the hmmm profiles are

outfolder=${alignmentsid_folder}/hmmsearch_${group}/


mkdir ${outfolder}

for i in $(ls ${alignmentsid_folder}/*.fa_oneline)
do
    sed "s/-//g" ${i} > ${outfolder}/${i}_nogap;
done #removing the gaps from the alignments to perform the hmm search


cat ${hmm_folder}/*.hmm  > ${outfolder}/hmm_profiles_${group}.hmm

hmmpress ${outfolder}/hmm_profiles_${group}.hmm

for j in ${outfolder}/*_nogap
do
    b=$(basename ${j})
    hmmscan -E 1e-10 --tblout ${outfolder}/${b}_resfile.txt ${outfolder}/hmm_profiles_${group}.hmm ${j} > /dev/null
done #performs the hmmmsearch



#cd resfiles

#min_sizes=$(stat -c "%s" *_resfile.txt | sort -n -r | tail -n 1)

#find . -name "*_resfile.txt" -size ${min_sizes}c -delete #deletes files with no matches based on the minimal size of iles



#grep "^POG09.*\s" -m 1 * | awk '$3 == 0' > best_matches.txt #select best matches for each file and delete the one that have low scores (ie do not match the first entry which should be 0

#awk '{print $1,$5}' best_matches.txt | sort -g -k2 | sed "s/:/ /" | awk '!seen[$2]++' > best_unique_matches.txt
#selects name of match and evalue, sorts the evalue, formats the column, select only the first best match of each profile
