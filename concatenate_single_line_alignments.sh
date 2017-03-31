target_folder=$1

paste -d'\0' ${target_folder}/*_oneline | sed "s/>\([0-9]\+\)>.*/>\1/"  | sed "/^\s*$/d" > concatenated_genes.fa

