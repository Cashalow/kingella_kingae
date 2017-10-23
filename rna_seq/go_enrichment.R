library(goseq)
library(stringr)
library(RMySQL)

calculate_goseq <- function(genes, gene_lengths, categories){
    vector <- as.integer(names(gene_lengths)%in%genes)
    names(vector) <- names(gene_lengths)
    pwf <- nullp(vector, bias.data=gene_lengths)
    res <- goseq(pwf, gene2cat=categories)
    return(res)   
}

extract_significant_cat <- function(res){
    return(res[res[,"over_represented_pvalue"]<0.05,])
}

m <- dbDriver("MySQL")
custom_tables <- dbConnect(m, group="custom_tables")
COG <- dbConnect(m, group="COG")
enzyme <- dbConnect(m, group="enzyme")

query_gos <- "SELECT * from uniprot_go_terms_chlamydia_04_16 LEFT JOIN locus2seqfeature_id_chlamydia_04_16 ON uniprot_go_terms_chlamydia_04_16.seqfeature_id = locus2seqfeature_id_chlamydia_04_16.seqfeature_id WHERE locus2seqfeature_id_chlamydia_04_16.locus_tag LIKE \"WCW_RS%\";"
query_cogs <- "SELECT DISTINCT description, functon, locus_tag2gi_hit_chlamydia_04_16.locus_tag, cog_names_2014.COG_id, cog_names_2014.name FROM locus_tag2gi_hit_chlamydia_04_16 LEFT JOIN (cog_names_2014, code2category) ON (cog_names_2014.COG_id = locus_tag2gi_hit_chlamydia_04_16.COG_id AND code2category.code = cog_names_2014.functon) WHERE locus_tag2gi_hit_chlamydia_04_16.locus_tag LIKE \"WCW_RS%\";"
query_keggs <- "SELECT * FROM pathway2ko LEFT JOIN (locus2ko_chlamydia_04_16, kegg_pathway) ON (locus2ko_chlamydia_04_16.ko_id = pathway2ko.ko_id and kegg_pathway.pathway_id = pathway2ko.pathway_id) where locus_tag like \"WCW_RS%\";"



all_gos <- dbGetQuery(custom_tables, query_gos)
gos_wcw <- data.frame(cbind(all_gos[,"locus_tag"], all_gos["go_description"]))

all_cog <- dbGetQuery(COG, query_cogs)
cog_wcw <- data.frame(cbind(all_cog[,"locus_tag"], all_cog["description"]))

all_keggs <- dbGetQuery(enzyme, query_keggs)
kegg_wcw <- data.frame(cbind(all_keggs[,"locus_tag"], all_keggs["description"]))

#read file with lengths

all_data <- read.csv("/home/sacha/Documents/rna_seq/2017-04-27_LGR-6-17_Wchondrophila-DeSeq2.txt", stringsAsFactor=F, sep="\t")
rownames(all_data) <- str_extract(all_data[,"Description"], "WCW_RS[0-9]+")

only24144 <- read.csv("/home/sacha/Documents/rna_seq/2017-04-28_LGR-12-17_Wchondrophila-DeSeq2.txt", stringsAsFactor=F, sep="\t")
rownames(only24144) <- str_extract(only24144[,"Description"], "WCW_RS[0-9]+")

h24 <- rownames(only24144[only24144[,"padj_144h.24h"]<0.05 & only24144[,"log2FoldChange_144h.24h"]<0,])
h144 <- rownames(only24144[only24144[,"padj_144h.24h"]<0.05 & only24144[,"log2FoldChange_144h.24h"]>0,])

h24p001 <- rownames(only24144[only24144[,"padj_144h.24h"]<0.01 & only24144[,"log2FoldChange_144h.24h"]<0,])
h144p001 <- rownames(only24144[only24144[,"padj_144h.24h"]<0.01 & only24144[,"log2FoldChange_144h.24h"]>0,])



values <- matrix(as.numeric(unlist(all_data[, 25:33])), ,9)
colnames(values) <- colnames(all_data[, 25:33])
rownames(values) <- rownames(all_data)

def <- str_extract(all_data[all_data[,"padj_6h.6h.def"]<0.01 & all_data[,"log2FoldChange_6h.6h.def"]<0,"Description"], "WCW_RS[0-9]+")


h6p001 <- str_extract(all_data[all_data[,"padj_24h.6h"]<0.01 & all_data[,"log2FoldChange_24h.6h"]<0,"Description"], "WCW_RS[0-9]+")
h6p001 <- h6p001[!is.na(h6p001)]

h6vs144p001 <- str_extract(all_data[all_data[,"padj_144h.6h"]<0.01 & all_data[,"log2FoldChange_144h.6h"]<0,"Description"], "WCW_RS[0-9]+")
h6vs144p001 <- h6vs144p001[!is.na(h6vs144p001)]


gene_lengths <- all_data[,5] - all_data[,4]
names(gene_lengths) <- str_extract(all_data[,"Description"], "WCW_RS[0-9]+")


#read corres between locus tag and seqfeature id

read_horn <- read.csv("/home/sacha/Documents/rna_seq/horn_res.csv", sep="\t", header=TRUE, stringsAsFactors=F)
res_horn<- matrix(as.numeric(unlist(read_horn[,c("X2.hpi_1", "X2.hpi_2", "X2.hpi_3", "X48.hpi_1", "X48.hpi_2" ,"X48.hpi_3", "X96.hpi_1" , "X96.hpi_2", "X96.hpi_3", "extracellular_1" , "extracellular_2"  , "extracellular_3" )])), ,12)
res_horn <- res_horn[complete.cases(res_horn),]
scaled_res <-  scale(res_horn, center=TRUE, scale=FALSE)
means <- c()
for (i in 0:3){
    means<-cbind(means, rowMeans(scaled_res[,(i*3+1):((i+1)*3)]))    
}

means2 <- c()
for (i in 0:3){
    means2<-cbind(means2, rowMeans(res_horn[,(i*3+1):((i+1)*3)]))    
}
scaled_res2 <- scale(means2, center=TRUE, scale=FALSE)

allgo <- read.csv("/home/sacha/Documents/rna_seq/allgo.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("seqID", "GO", "description", "locus_tag", "seqID2", "taxID"))
gos_wcw <- data.frame(cbind(allgo[,"locus_tag"], allgo["description"]))

allkegg <- read.csv("/home/sacha/Documents/rna_seq/kegg.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("pathway_id", "keggID", "taxon_id", "locus_tag", "orthogroup", "keggID2","pathway_id2", "pathway_name", "short", "category", "general"))
kegg_wcw <- data.frame(cbind(allkegg[,"locus_tag"], allkegg["general"]))

allcog <- read.csv("/home/sacha/Documents/rna_seq/cog.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("category","locus_tag", "COG_ID", "description"))
cog_wcw <- data.frame(cbind(allcog[,"locus_tag"], allcog["category"]))

gosh6vs24 <- calculate_goseq(h6p001, gene_lengths, gos_wcw)
gosh6vs144 <- calculate_goseq(h6vs144p001, gene_lengths, gos_wcw)

gosh24 <- calculate_goseq(h24, gene_lengths, gos_wcw)
gosh144 <- calculate_goseq(h144, gene_lengths, gos_wcw)

gosh24p001 <- calculate_goseq(h24p001, gene_lengths, gos_wcw)
gosh144p001 <- calculate_goseq(h144p001, gene_lengths, gos_wcw)

for (i in list(gosh24, gosh24p001, gosh144, gosh144p001)){
    print(extract_significant_cat(i)[,c("category","over_represented_pvalue")])
    }





keggh6vs24 <- calculate_goseq(h6p001, gene_lengths, kegg_wcw)
keggh6vs144 <- calculate_goseq(h6vs144p001, gene_lengths, kegg_wcw)

keggh24 <- calculate_goseq(h24, gene_lengths, kegg_wcw)
keggh144 <- calculate_goseq(h144, gene_lengths, kegg_wcw)
keggh24p001 <- calculate_goseq(h24p001, gene_lengths, kegg_wcw)
keggh144p001 <- calculate_goseq(h144p001, gene_lengths, kegg_wcw)


for (i in list(keggh24, keggh24p001, keggh144, keggh144p001)){
    print(extract_significant_cat(i)[,c("category","over_represented_pvalue")])
    }




cogh6 <- calculate_goseq(h6p001, gene_lengths, cog_wcw)

cogh24 <- calculate_goseq(h24, gene_lengths, cog_wcw)
cogh144 <- calculate_goseq(h144, gene_lengths, cog_wcw)
cogh24p001 <- calculate_goseq(h24p001, gene_lengths, cog_wcw)
cogh144p001 <- calculate_goseq(h144p001, gene_lengths, cog_wcw)

for (i in list(cogh24, cogh24p001, cogh144, cogh144p001)){
    print(extract_significant_cat(i)[,c("category","over_represented_pvalue")])
    }
