orthoselect=$1

#this is thefile with the names of the orthogroups to extract

ids=/home/sacha/Documents/Kingella/comparative_genomics/orthofinder/kingella_neisseria/Results_Mar17/WorkingDirectory/SpeciesIDs.txt_cop 

orthofile=/home/sacha/Documents/Kingella/comparative_genomics/orthofinder/kingella_neisseria/Results_Mar17/Orthogroups.csv

destfolder=/home/sacha/Documents/Kingella/comparative_genomics/genetic_alignment/second/

annotationfolder=/home/sacha/Documents/Kingella/annotations/

mkdir $destfolder

for i in $( tail -n +2 ${ids} | cut -f1 -d' ' | sed "s/://")
do
    col=$(( i+2 ))
    line=$(( i+1 ))
    for j in $(cat ${orthoselect})
    do
	gi=$(grep $j ${orthofile} | cut -f ${col} | cut -d' ' -f1 )
    	echo ${i} ${gi} >> ${destfolder}/genelists/${j}.txt
	species=$(cut -f 2 -d' ' ${ids} | sed "${line}q;d" | sed 's/\.faa//')
	echo ">${i}" >> ${destfolder}/fastas/${j}.fa
	sed -n -e "/${gi}/,/^>/ p" ${annotationfolder}/${species}/${species}.ffn | sed -e '1d;$d' >> ${destfolder}/fastas/${j}.fa
    done
done

					
