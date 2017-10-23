library(goseq)
library(stringr)
library(RMySQL)
library(DESeq2)

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

#DEFINE GROUPS in ~/.my.conf to connect to the local MySQL DB (username, password, db name, host needed)

m <- dbDriver("MySQL")
custom_tables <- dbConnect(m, group="custom_tables")
COG <- dbConnect(m, group="COG")
enzyme <- dbConnect(m, group="enzyme")

query_gos <- "SELECT * from uniprot_go_terms_chlamydia_04_16 LEFT JOIN locus2seqfeature_id_chlamydia_04_16 ON uniprot_go_terms_chlamydia_04_16.seqfeature_id = locus2seqfeature_id_chlamydia_04_16.seqfeature_id WHERE locus2seqfeature_id_chlamydia_04_16.locus_tag LIKE \"WCW_RS%\";"

query_cogs <- "SELECT DISTINCT description, functon, locus_tag2gi_hit_chlamydia_04_16.locus_tag, cog_names_2014.COG_id, cog_names_2014.name FROM locus_tag2gi_hit_chlamydia_04_16 LEFT JOIN (cog_names_2014, code2category) ON (cog_names_2014.COG_id = locus_tag2gi_hit_chlamydia_04_16.COG_id AND code2category.code = cog_names_2014.functon) WHERE locus_tag2gi_hit_chlamydia_04_16.locus_tag LIKE \"WCW_RS%\";"

query_keggs <- "SELECT * FROM pathway2ko LEFT JOIN (locus2ko_chlamydia_04_16, kegg_pathway) ON (locus2ko_chlamydia_04_16.ko_id = pathway2ko.ko_id AND kegg_pathway.pathway_id = pathway2ko.pathway_id) WHERE locus_tag LIKE \"WCW_RS%\";"



all_gos <- dbGetQuery(custom_tables, query_gos)
gos_wcw <- data.frame(cbind(all_gos[,"locus_tag"], all_gos["go_description"]))

all_cogs <- dbGetQuery(COG, query_cogs)
cog_wcw <- data.frame(cbind(all_cogs[,"locus_tag"], all_cogs["description"]))

all_keggs <- dbGetQuery(enzyme, query_keggs)
kegg_wcw <- data.frame(cbind(all_keggs[,"locus_tag"], all_keggs["description"]))

#read file with lengths

all_data <- read.csv("/home/sacha/Documents/rna_seq/2017-04-27_LGR-6-17_Wchondrophila-DeSeq2.txt", stringsAsFactor=F, sep="\t")
rownames(all_data) <- str_extract(all_data[,"Description"], "WCW_RS[0-9]+")

only24144 <- read.csv("/home/sacha/Documents/rna_seq/2017-04-28_LGR-12-17_Wchondrophila-DeSeq2.txt", stringsAsFactor=F, sep="\t")
rownames(only24144) <- str_extract(only24144[,"Description"], "WCW_RS[0-9]+")

genes <- read.csv("/home/sacha/Documents/rna_seq/genes.csv", stringsAsFactor=F, sep="\t", row.names=1)

#get enriched genes
h24p005 <- rownames(only24144[only24144[,"padj_144h.24h"]<0.05 & only24144[,"log2FoldChange_144h.24h"]<0,])
h144p005 <- rownames(only24144[only24144[,"padj_144h.24h"]<0.05 & only24144[,"log2FoldChange_144h.24h"]>0,])

h24p001 <- rownames(only24144[only24144[,"padj_144h.24h"]<0.01 & only24144[,"log2FoldChange_144h.24h"]<0,])
h144p001 <- rownames(only24144[only24144[,"padj_144h.24h"]<0.01 & only24144[,"log2FoldChange_144h.24h"]>0,])

h6vs24p001 <- rownames(all_data[all_data[,"padj_24h.6h"]<0.01 & all_data[,"log2FoldChange_24h.6h"]>0 & !is.na(all_data[,"padj_24h.6h"]),])

h6vs144p001 <- rownames(all_data[all_data[,"padj_144h.6h"]<0.01 & all_data[,"log2FoldChange_144h.6h"]>0 & !is.na(all_data[,"padj_144h.6h"]),])

h6vs6def <- rownames(all_data[all_data[,"padj_6h.6h.def"]<0.05 & all_data[,"log2FoldChange_6h.6h.def"]<0 & !is.na(all_data[,"padj_6h.6h.def"]),])

h24vs6def <- rownames(all_data[all_data[,"padj_24h.6h.def"]<0.05 & all_data[,"log2FoldChange_24h.6h.def"]<0 & !is.na(all_data[,"padj_24h.6h.def"]),])

#redo DeSeq2 analysis

raw_counts <- all_data[,10:21]

groups <- data.frame(timepoint=c(rep("6h",2), rep("24h", 3), rep("144h", 2), "24h", rep("6h", 3), "6h"), treatment=c(rep("def",2), rep("none",9), "def"))
rownames(groups) <- colnames(all_data)[10:21]
raw_counts <- raw_counts[,rownames(groups)]

lib_notreat <- c("LGR.12","LGR.13","LGR.14","LGR.15","LGR.16","LGR.17","LGR.6","LGR.7","LGR.8")
lib_h14424 <- c("LGR.12","LGR.13","LGR.14","LGR.15","LGR.16","LGR.17")
lib_h6 <- c("LGR.6","LGR.7","LGR.8", "LGR.10", "LGR.11", "LGR.9")

groups_h6 <- data.frame(treatment=groups[lib_h6,"treatment"])
colnames(groups_h6) <- "treatment"

dds<-DESeqDataSetFromMatrix(countData=raw_counts[,lib_notreat] , colData=data.frame(timepoint=groups[lib_notreat,"timepoint"]), design=~timepoint)
results <- DESeq(dds)

h6 <- results(results, name="timepoint6h")
h24 <- results(results, name="timepoint24h")
h144 <- results(results, name="timepoint144h")
    
resh6 <- rownames(h6[h6[,"padj"] < 0.05 & h6[,"log2FoldChange"]>0 & !is.na(h6[,"padj"]),])
resh24 <- rownames(h24[h24[,"padj"] < 0.05 & h24[,"log2FoldChange"]>0 & !is.na(h24[,"padj"]),])
resh144 <- rownames(h144[h144[,"padj"] < 0.05 & h144[,"log2FoldChange"]>0,])



extract_significant_cat(calculate_goseq(resh6, gene_lengths, gos_wcw))
extract_significant_cat(calculate_goseq(resh24, gene_lengths, gos_wcw))
extract_significant_cat(calculate_goseq(resh144, gene_lengths, gos_wcw))



dds<-DESeqDataSetFromMatrix(countData=raw_counts[,lib_h14424] , colData=groups[lib_h14424,], design=~timepoint)
results <- DESeq(dds)

h24 <- results(results, name="timepoint24h")
h144 <- results(results, name="timepoint144h")
    
resh144 <- rownames(h144[h144[,"padj"] < 0.05 & h144[,"log2FoldChange"]>0,])
resh24 <- rownames(h24[h24[,"padj"] < 0.01 & h24[,"log2FoldChange"]>0 & !is.na(h24[,"padj"]),])

extract_significant_cat(calculate_goseq(resh144, gene_lengths, gos_wcw))
extract_significant_cat(calculate_goseq(resh24, gene_lengths, gos_wcw))


dds<-DESeqDataSetFromMatrix(countData=raw_counts[,lib_h6] , colData=data.frame(treatment=groups[lib_h6,"treatment"]), design=~treatment)
results <- DESeq(dds)

def <- results(results, name="treatmentdef")
    
resdef <- rownames(def[def[,"pvalue"] < 0.05 & def[,"log2FoldChange"]>0 & !is.na(def[,"pvalue"]),])

extract_significant_cat(calculate_goseq(resdef, gene_lengths, gos_wcw))

#get gene lengths from result file (should be done from DB)
gene_lengths <- all_data[,5] - all_data[,4]
names(gene_lengths) <- rownames(all_data)

#get enriched categories

#GOs
gosh6vs24 <- calculate_goseq(h6vs24p001, gene_lengths, gos_wcw)
gosh6vs144 <- calculate_goseq(h6vs144p001, gene_lengths, gos_wcw)

gosh24 <- calculate_goseq(h24p005, gene_lengths, gos_wcw)
gosh144 <- calculate_goseq(h144p005, gene_lengths, gos_wcw)

gosh24p001 <- calculate_goseq(h24p001, gene_lengths, gos_wcw)
gosh144p001 <- calculate_goseq(h144p001, gene_lengths, gos_wcw)

for (i in list(gosh6vs24, gosh6vs144, gosh24, gosh24p001, gosh144, gosh144p001)){
    print(extract_significant_cat(i)[,c("category","over_represented_pvalue")])
}


#KEGG
keggh6vs24 <- calculate_goseq(h6p001, gene_lengths, kegg_wcw)
keggh24 <- calculate_goseq(h24p005, gene_lengths, kegg_wcw)
keggh144 <- calculate_goseq(h144p005, gene_lengths, kegg_wcw)

keggh6vs144 <- calculate_goseq(h6vs144p001, gene_lengths, kegg_wcw)
keggh24p001 <- calculate_goseq(h24p001, gene_lengths, kegg_wcw)
keggh144p001 <- calculate_goseq(h144p001, gene_lengths, kegg_wcw)


for (i in list(keggh6vs24, keggh6vs144, keggh24, keggh24p001, keggh144, keggh144p001)){
    print(extract_significant_cat(i)[,c("category","over_represented_pvalue")])
}



#COGs

cogh6 <- calculate_goseq(h6p001, gene_lengths, cog_wcw)
cogh24 <- calculate_goseq(h24p005, gene_lengths, cog_wcw)
cogh144 <- calculate_goseq(h144p005, gene_lengths, cog_wcw)

cogh24p001 <- calculate_goseq(h24p001, gene_lengths, cog_wcw)
cogh144p001 <- calculate_goseq(h144p001, gene_lengths, cog_wcw)

for (i in list(cogh24, cogh24p001, cogh144, cogh144p001)){
    print(extract_significant_cat(i)[,c("category","over_represented_pvalue")])
}



#work on expression values
all_exp_values <- matrix(as.numeric(unlist(all_data[, 25:33])), ,9)
colnames(all_exp_values) <- colnames(all_data[, 25:33])
rownames(all_exp_values) <- rownames(all_data)

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
