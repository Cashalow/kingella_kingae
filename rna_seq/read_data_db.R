library(RMySQL)

#DEFINE GROUPS in ~/.my.conf to connect to the local MySQL DB (username, password, db name, host needed)

m <- dbDriver("MySQL")
custom_tables <- dbConnect(m, group="custom_tables")
COG <- dbConnect(m, group="COG")
enzyme <- dbConnect(m, group="enzyme")

query_gos <- "SELECT * from uniprot_go_terms_chlamydia_04_16 LEFT JOIN locus2seqfeature_id_chlamydia_04_16 ON uniprot_go_terms_chlamydia_04_16.seqfeature_id = locus2seqfeature_id_chlamydia_04_16.seqfeature_id WHERE locus2seqfeature_id_chlamydia_04_16.locus_tag LIKE \"WCW_RS%\";"

query_cogs <- "SELECT DISTINCT description, functon, locus_tag2gi_hit_chlamydia_04_16.locus_tag, cog_names_2014.COG_id, cog_names_2014.name FROM locus_tag2gi_hit_chlamydia_04_16 LEFT JOIN (cog_names_2014, code2category) ON (cog_names_2014.COG_id = locus_tag2gi_hit_chlamydia_04_16.COG_id AND code2category.code = cog_names_2014.functon) WHERE locus_tag2gi_hit_chlamydia_04_16.locus_tag LIKE \"WCW_RS%\";"

query_keggs <- "SELECT * FROM pathway2ko LEFT JOIN (locus2ko_chlamydia_04_16, kegg_pathway) ON (locus2ko_chlamydia_04_16.ko_id = pathway2ko.ko_id AND kegg_pathway.pathway_id = pathway2ko.pathway_id) WHERE locus_tag LIKE \"WCW_RS%\";"

query_uniprot_annot <- "SELECT * from uniprot_annotation_chlamydia_04_16 LEFT JOIN locus2seqfeature_id_chlamydia_04_16 ON uniprot_annotation_chlamydia_04_16.seqfeature_id=locus2seqfeature_id_chlamydia_04_16.seqfeature_id WHERE locus2seqfeature_id_chlamydia_04_16.locus_tag like \"WCW_RS%\";"



all_gos <- dbGetQuery(custom_tables, query_gos)
gos_wcw <- data.frame(cbind(all_gos[,"locus_tag"], all_gos["go_description"]))

all_cogs <- dbGetQuery(COG, query_cogs)
cog_wcw <- data.frame(cbind(all_cogs[,"locus_tag"], all_cogs["description"]))

all_keggs <- dbGetQuery(enzyme, query_keggs)
kegg_wcw <- data.frame(cbind(all_keggs[,"locus_tag"], all_keggs["description"]))

genes <- dbGetQuery(custom_tables, query_uniprot_annot)
gene_names <- data.frame(names=genes[,"gene"])
rownames(gene_names) <- genes[,"locus_tag"]

