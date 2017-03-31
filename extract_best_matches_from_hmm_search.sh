hmm_res_folder=$1
preffix=$2 #this is the preffix of hmm profiles to search for

grep "^${preffix}" -L ${hmm_res_folder}/* | xargs -I % rm % #deletes files with no matches


grep "^${preffix}.*\s" -m 1 * | awk '$3 == 0' > ${hmm_res_folder}/best_matches.txt #select best matches for each file and delete the one that have low scores (ie do not match the first entry which should be 0

awk '{print $1,$5}' ${hmm_res_folder}/best_matches.txt | sort -g -k2 | sed "s/:/ /" | awk '!seen[$2]++' | awk '{print $1}' | sed "s/_oneline.*//" > ${hmm_res_folder}/best_unique_matches.txt
#selects name of match and evalue, sorts the evalue, formats the column, select only the first best match of each profile

