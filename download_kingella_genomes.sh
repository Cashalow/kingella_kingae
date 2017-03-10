kingella_genome_list_adress="https://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=download&orgn=%22Kingella%20kingae%22[orgn]&status=50|30|20&report=proks&group=--%20All%20Prokaryotes%20--&subgroup=--%20All%20Prokaryotes%20--&format="

wget $kingella_genome_list_adress -O kingella_genome_list.txt

tail -n +2 kingella_genome_list.txt | cut -f20 | sed -s "s/\(.*\)\/\(.*\)$/\1\/\2\/\2\_genomic.fna.gz/" | sed -s "s/\r//g" > kingella_genome_adresses.txt

mkdir raw_genome_sequences
cd raw_genome_sequences
wget -i ../kingella_genome_adresses.txt
gzip -d * 


tail -n +2 kingella_genome_list.txt | cut -f20 | sed -s "s/\(.*\)\/\(.*\)$/\1\/\2\/\2\_genomic.gff.gz/" | sed -s "s/\r//g" > kingella_annotations_adresses.txt
