library(goseq)
library(stringr)

#read file with lengths

all_data <- read.csv("2017-04-28_LGR-12-17_Wchondrophila-DeSeq2.txt", stringsAsFactor=F, sep="\t")


h24 <- str_extract(all_data[all_data[,"padj_144h.24h"]<0.05 & all_data[,"log2FoldChange_144h.24h"]<0,"Description"], "WCW_RS[0-9]+")

h144 <- str_extract(all_data[all_data[,"padj_144h.24h"]<0.05 & all_data[,"log2FoldChange_144h.24h"]>0,"Description"], "WCW_RS[0-9]+")


gene_lengths <- all_data[,5] - all_data[,4]
names(gene_lengths) <- str_extract(all_data[,"Description"], "WCW_RS[0-9]+")

#read corres between locus tag and seqfeature id


allgo <- read.csv("allgo.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("seqID", "GO", "description", "locus_tag", "seqID2", "taxID"))

gos_wcw <- data.frame(cbind(allgo[,"locus_tag"], allgo["description"]))

vectorh24 <- as.integer(names(gene_lengths)%in%h24)
names(vectorh24) <- names(gene_lengths)
pwfh24 <- nullp(vectorh24, bias.data=gene_lengths)
resh24 <- goseq(pwfh24, gene2cat=gos_wcw)

vectorh144 <- as.integer(names(gene_lengths)%in%h144)
names(vectorh144)<-names(gene_lengths)
pwfh144 <- nullp(vectorh144, bias.data=gene_lengths)
resh144 <- goseq(pwfh144, gene2cat=gos_wcw)

allkegg <- read.csv("kegg.csv", sep="\t", header=FALSE, stringsAsFactors=F, col.names=c("keggID", "locus_tag", "KO", "category", "general", "description"))

kegg_wcw <- data.frame(cbind(allkegg[,"locus_tag"], allkegg["general"]))


kegg24 <- goseq(pwfh24, gene2cat=kegg_wcw)
kegg144 <- goseq(pwfh144, gene2cat=kegg_wcw)
