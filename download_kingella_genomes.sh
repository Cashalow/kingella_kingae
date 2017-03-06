kingella_genome_list_adress="https://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=download&orgn=%22Kingella%20kingae%22[orgn]&status=50|30|20&report=proks&group=--%20All%20Prokaryotes%20--&subgroup=--%20All%20Prokaryotes%20--&format="

wget $kingella_genome_list_adress -o kingella_genome_list.txt
cut -f19 kigella_genome_list.txt | sed -s "s/\(.*\)\/\(.*\)$/\1\/\2\/\2_genomic.fna.gz/" > kingella_genome_adresses.txt


