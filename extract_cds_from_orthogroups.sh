orthoselect=$1

#this is the file with the names of the orthogroups to extract

ids=$2
#file with the species ID/names

orthofile=$3
#csv file of the orthogroups 

destfolder=$4

annotationfolder=$5
#where the cds are

mkdir -p $destfolder/genelists/
mkdir -p $destfolder/fastas/

for i in $( cat ${ids} | cut -f1 -d' ' | sed "s/://")
do
    echo ${i}
    col=$(( i+2 ))
    line=$(( i+1 ))
    for j in $(cat ${orthoselect} | sed "s/.fa//")
    do
	gi=$(grep $j ${orthofile} | cut -f ${col} | cut -d' ' -f1 )
    	echo ${i} ${gi} >> ${destfolder}/genelists/${j}.txt
	species=$(cut -f 2 -d' ' ${ids} | sed "${line}q;d" | sed 's/\.faa//')
	echo ">${i}" >> ${destfolder}/fastas/${j}.fa
	sed -n -e "/${gi}/,/^>/ p" ${annotationfolder}/${species}/${species}.ffn | sed -e '1d;$d' >> ${destfolder}/fastas/${j}.fa
    done
done

					
