file_with_ogs=$1
alignmentsid_folder=$2 #folder with the alignments and numbers as sole ids for species

output_folder=$( dirname ${file_with_ogs} )

for og in $(cat ${file_with_ogs})
do
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${alignmentsid_folder}/${og} > ${output_folder}/${og}_oneline;
done

sed -i "s/>\([0-9]\+\)_.*/>\1/" *_oneline
sed -i "/^\s*$/d" *_oneline


