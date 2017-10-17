library(goseq)
library(stringr)

calculate_goseq <- function(genes, gene_lengths, categories){
    vector <- as.integer(names(gene_lengths)%in%genes)
    names(vector) <- names(gene_lengths)
    pwf <- nullp(vector, bias.data=gene_lengths)
    res <- goseq(pwf, gene2cat=categories)
    return(res)   
}

#read file with lengths

all_data <- read.csv("/home/sacha/Documents/rna_seq/2017-04-27_LGR-6-17_Wchondrophila-DeSeq2.txt", stringsAsFactor=F, sep="\t")
    
h24 <- str_extract(all_data[all_data[,"padj_144h.24h"]<0.05 & all_data[,"log2FoldChange_144h.24h"]<0,"Description"], "WCW_RS[0-9]+")

h144 <- str_extract(all_data[all_data[,"padj_144h.24h"]<0.05 & all_data[,"log2FoldChange_144h.24h"]>0,"Description"], "WCW_RS[0-9]+")

h6 <- str_extract(all_data[all_data[,"padj_24h.6h"]<0.05 & all_data[,"log2FoldChange_24h.6h"]<0,"Description"], "WCW_RS[0-9]+")
h6 <- h6[!is.na(h6)]

h6vs144 <- str_extract(all_data[all_data[,"padj_144h.6h"]<0.05 & all_data[,"log2FoldChange_144h.6h"]<0,"Description"], "WCW_RS[0-9]+")
h6vs144 <- h6vs144[!is.na(h6vs144)]


gene_lengths <- all_data[,5] - all_data[,4]
names(gene_lengths) <- str_extract(all_data[,"Description"], "WCW_RS[0-9]+")


#read corres between locus tag and seqfeature id



allgo <- read.csv("/home/sacha/Documents/rna_seq/allgo.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("seqID", "GO", "description", "locus_tag", "seqID2", "taxID"))
gos_wcw <- data.frame(cbind(allgo[,"locus_tag"], allgo["description"]))

allkegg <- read.csv("/home/sacha/Documents/rna_seq/kegg.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("keggID", "locus_tag", "KO", "category", "general", "description"))
kegg_wcw <- data.frame(cbind(allkegg[,"locus_tag"], allkegg["general"]))

allcog <- read.csv("/home/sacha/Documents/rna_seq/cog.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("category","locus_tag", "COG_ID", "description"))
cog_wcw <- data.frame(cbind(allcog[,"locus_tag"], allcog["category"]))

gosh6 <- calculate_goseq(h6, gene_lengths, gos_wcw)
gosh24 <- calculate_goseq(h24, gene_lengths, gos_wcw)
gosh144 <- calculate_goseq(h144, gene_lengths, gos_wcw)

keggh6 <- calculate_goseq(h6, gene_lengths, kegg_wcw)
keggh24 <- calculate_goseq(h24, gene_lengths, kegg_wcw)
keggh144 <- calculate_goseq(h144, gene_lengths, kegg_wcw)

cogh6 <- calculate_goseq(h6, gene_lengths, cog_wcw)
cogh24 <- calculate_goseq(h24, gene_lengths, cog_wcw)
cogh144 <- calculate_goseq(h144, gene_lengths, cog_wcw)
