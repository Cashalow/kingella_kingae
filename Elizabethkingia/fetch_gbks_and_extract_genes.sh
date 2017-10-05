#/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
esearch -db assembly -query "Elizabethkingia" | efetch -db assembly -format docsum | xtract -pattern DocumentSummary -if RefSeq -element SpeciesName FtpPath_RefSeq > refseq_ftp.txt
cut -f 2 refseq_ftp.txt  | sed "s/\(\/GCF_.*\)/\1\1_genomic.gbff.gz/" > refseq_genbanks.txt
mkdir -p gbks && cd gbks
while read p;
do
    b=$(basename ${p})
    if [ ! -f ${b%.gz}  ]
    then
	wget ${p}
	gzip -f -d ${b}
    fi
done < ../refseq_genbanks.txt
mkdir -p fasta && cd fasta 
rm -f *.fasta && for i in $(ls ../*.gbff); do python ${DIR}/extract_gbk.py ${i}; done
