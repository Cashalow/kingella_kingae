list_of_files <- c("res_yverdon", "res_oralis", "res_denitrificans", "res_potus")

results <- list()

for (i in list_of_files){
    results[[i]] <- read.table(i, header=FALSE, stringsAsFactors=FALSE, sep="\t")
    colnames(results[[i]]) <- c("gene", paste(i, c("id", "cov", "evalue")))
}
    
Reduce(function(x, y) merge(x, y, by="gene", all=T), lapply(results, function(x) x[x[,4] < evalue,]))[,c(1, 0:4*3+2)]
