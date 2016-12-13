


rm(list=ls())

library(piano)
library(xlsx)
library(WriteXLS)
library(GSA)


aa <- read.xlsx("../Data/genePathways-allKEGG&MouseCyc.xlsx",1)

for (sourceName in c("KEGG", "MouseCyc")) {
 aa <- aa[aa$source==sourceName,]
  
 aa$pathway <- as.character(aa$pathway)
 aa$hgnc_symbol_ids <- toupper(as.character(aa$hgnc_symbol_ids))
 gTogs <- data.frame()
 for (i in 1:dim(aa)[1]) {
   x <- unlist(strsplit(aa$hgnc_symbol_ids[i], ","))
   gTogs <- rbind(gTogs,(cbind(x, GSetName=paste(aa$pathway[i], aa$source[i], sep="%"))))
 }
 length(unique(gTogs$GSetName))
 a <- loadGSC(gTogs)

 d <- read.table("../Data/gene_exp.diff")
 d2 <- subset(d, d$V5=='K' & d$V6 =='L' & d$V11 !="-nan") ## NOTE --> 23341
 ## filter KRT
 indices <- gdata::startsWith(as.character(d2$V2), "Krt")
 if (sum(indices)==0) stop("error!")
 d2 <- d2[!indices,]
 universe <- d2$V1

 universe <- toupper(universe)
 ## subset the significant genes
 sig <- d2[d2$V14=="yes",] ## 1714
 sig$V10 <- as.numeric(as.character(sig$V10))

 ## OVER-REP analysis the genes UP-regulated in K
 sigUpinK <- sig[sig$V10< -1.25,] ## 915
 genes <- sigUpinK$V1  ## 
 genes <- toupper(genes)
 significant <- rep(1, length(universe))
 significant[which(universe%in%genes)] <- 0
 res <- runGSAhyper(genes=universe, significant,  gsc=a, adjMethod="fdr")
 result <- res$resTab[order(res$resTab[, "p-value"]),]
 Sys.setlocale('LC_ALL','C') 
 write.xlsx2(result, file= paste("../Ouput/Fig3/KL-allSig_", sourceName, "-UP.xlsx", sep=""))

 ## OVER-REP analysis the genes DOWN-regulated in K
 sigDninK <- sig[sig$V10> 1.25,] ## 799
 genes <- sigDninK$V1
 genes <- toupper(genes)
 significant <- rep(1, length(universe))
 significant[which(universe%in%genes)] <- 0
 res <- runGSAhyper(genes=universe, significant,  gsc=a, adjMethod="fdr")
 result <- res$resTab[order(res$resTab[, "p-value"]),]
 Sys.setlocale('LC_ALL','C') 
 write.xlsx2(result, file= paste("../Ouput/Fig3/KL-allSig_", sourceName, "-DN.xlsx", sep=""))
}


