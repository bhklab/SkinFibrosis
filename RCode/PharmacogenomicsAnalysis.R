

## Core of the computational pharmacogenomis analysis pipeline has been impelemented in this script 
## the script runs for 6 different query sizes: from 150 to 400, intervals of 50
## the script generates material for Fig4a-PharmacoGenomics.R

rm(list=ls())

library(PharmacoGx)
library(Biobase)
library(parallel)
library(xlsx)

## this object contains the (most significant) diferentially expressed genes at the last time point between
##the following experimental groups:  (irradiated+) adsc-treated and irradiated 
load("../Data/useBiomart.RData")

ind <- match(dSignToHuman$hgnc_symbol, toupper(dSign$V1))
dSignToHuman <- cbind(dSignToHuman, ind)
dSignToHuman <- dSignToHuman[order(dSignToHuman$ind),]


for (n in seq(from=75,to=200, by=25)) {
   #### HEAD PART of the query
   qh <- head(dSignToHuman,n)
   if (sum(dSign[qh$ind,"V11"]>0)!=n) stop("ERROR!")
   sigHEntrzh <- lapply(qh$entrezgene,function(x) paste("geneid.",x, sep=''))
   qSigh <- as.numeric(rep(-1,length(sigHEntrzh)))
   names(qSigh) <- sigHEntrzh

   ####TAIL PART of the query
   qt <- tail(dSignToHuman,n)
   if (sum(dSign[qt$ind,"V11"]<0)!=n) stop("ERROR!")
   sigHEntrzt <- lapply(qt$entrezgene,function(x) paste("geneid.",x, sep=''))
   qSigt <- as.numeric(rep(1,length(sigHEntrzt)))
   names(qSigt) <- sigHEntrzt

   qSig <- c(qSigh,qSigt)
   ## 2- query the CMAP (parallel code on the server)
   ncor <- detectCores()
   if (ncor>16) ncor <- 16
   cl <- makeCluster(ncor)
   load("../Data/CMAP_gene_signatures.RData")
   res <- parApply(CMAP.genePerturbations[,,c("tstat", "fdr")], 2, function(x, QSIG){
       return(PharmacoGx::connectivityScore(x=x,
                                       y=QSIG,
                                       method="gsea", nperm=1000))
   }, cl = cl, QSIG=qSig)

   stopCluster(cl)
   rownames(res) <- c("Connectivity", "P Value")
   res <- t(res)
   res <- res[order(res[,1], decreasing=T),]
   #save(list = ls(all.names = TRUE), file = paste("qSig",length(qSig),".RData", sep=""))
   write.xlsx2(res, file= paste("../Ouput/DrugList-qSig",length(qSig), ".xlsx", sep=""))
}


