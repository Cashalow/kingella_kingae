library(goseq)
library(stringr)

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

#read file with lengths

all_data <- read.csv("/home/sacha/Documents/rna_seq/2017-04-27_LGR-6-17_Wchondrophila-DeSeq2.txt", stringsAsFactor=F, sep="\t")

only24144 <- read.csv("/home/sacha/Documents/rna_seq/2017-04-28_LGR-12-17_Wchondrophila-DeSeq2.txt", stringsAsFactor=F, sep="\t")

h24 <- str_extract(only24144[only24144[,"padj_144h.24h"]<0.05 & only24144[,"log2FoldChange_144h.24h"]<0,"Description"], "WCW_RS[0-9]+")
h144 <- str_extract(only24144[only24144[,"padj_144h.24h"]<0.05 & only24144[,"log2FoldChange_144h.24h"]>0,"Description"], "WCW_RS[0-9]+")

h24p001 <- str_extract(only24144[only24144[,"padj_144h.24h"]<0.01 & only24144[,"log2FoldChange_144h.24h"]<0,"Description"], "WCW_RS[0-9]+")
h144p001 <- str_extract(only24144[only24144[,"padj_144h.24h"]<0.01 & only24144[,"log2FoldChange_144h.24h"]>0,"Description"], "WCW_RS[0-9]+")



h6p001 <- str_extract(all_data[all_data[,"padj_24h.6h"]<0.01 & all_data[,"log2FoldChange_24h.6h"]<0,"Description"], "WCW_RS[0-9]+")
h6p001 <- h6p001[!is.na(h6p001)]

h6vs144p001 <- str_extract(all_data[all_data[,"padj_144h.6h"]<0.01 & all_data[,"log2FoldChange_144h.6h"]<0,"Description"], "WCW_RS[0-9]+")
h6vs144p001 <- h6vs144p001[!is.na(h6vs144p001)]


gene_lengths <- all_data[,5] - all_data[,4]
names(gene_lengths) <- str_extract(all_data[,"Description"], "WCW_RS[0-9]+")


#read corres between locus tag and seqfeature id



allgo <- read.csv("/home/sacha/Documents/rna_seq/allgo.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("seqID", "GO", "description", "locus_tag", "seqID2", "taxID"))
gos_wcw <- data.frame(cbind(allgo[,"locus_tag"], allgo["description"]))

allkegg <- read.csv("/home/sacha/Documents/rna_seq/kegg.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("pathway_id", "keggID", "taxon_id", "locus_tag", "orthogroup", "keggID2","pathway_id2", "pathway_name", "short", "category", "general"))
kegg_wcw <- data.frame(cbind(allkegg[,"locus_tag"], allkegg["general"]))

allcog <- read.csv("/home/sacha/Documents/rna_seq/cog.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("category","locus_tag", "COG_ID", "description"))
cog_wcw <- data.frame(cbind(allcog[,"locus_tag"], allcog["category"]))

gosh6 <- calculate_goseq(h6, gene_lengths, gos_wcw)

gosh24 <- calculate_goseq(h24, gene_lengths, gos_wcw)
gosh144 <- calculate_goseq(h144, gene_lengths, gos_wcw)

gosh24p001 <- calculate_goseq(h24p001, gene_lengths, gos_wcw)
gosh144p001 <- calculate_goseq(h144p001, gene_lengths, gos_wcw)

for (i in list(gosh24, gosh24p001, gosh144, gosh144p001)){
    print(extract_significant_cat(i)[,c("category","over_represented_pvalue")])
    }


keggh6 <- calculate_goseq(h6, gene_lengths, kegg_wcw)

keggh24 <- calculate_goseq(h24, gene_lengths, kegg_wcw)
keggh144 <- calculate_goseq(h144, gene_lengths, kegg_wcw)
keggh24p001 <- calculate_goseq(h24p001, gene_lengths, kegg_wcw)
keggh144p001 <- calculate_goseq(h144p001, gene_lengths, kegg_wcw)


for (i in list(keggh24, keggh24p001, keggh144, keggh144p001)){
    print(extract_significant_cat(i)[,c("category","over_represented_pvalue")])
    }



cogh6 <- calculate_goseq(h6, gene_lengths, cog_wcw)

cogh24 <- calculate_goseq(h24, gene_lengths, cog_wcw)
cogh144 <- calculate_goseq(h144, gene_lengths, cog_wcw)
cogh24p001 <- calculate_goseq(h24p001, gene_lengths, cog_wcw)
cogh144p001 <- calculate_goseq(h144p001, gene_lengths, cog_wcw)

for (i in list(cogh24, cogh24p001, cogh144, cogh144p001)){
    print(extract_significant_cat(i)[,c("category","over_represented_pvalue")])
    }
