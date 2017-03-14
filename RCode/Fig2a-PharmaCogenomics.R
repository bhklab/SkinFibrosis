
rm(list=ls())

## this script performs analysis required to generate fig2A
source(./RCode/querySig.R)
source(./RCode/plotDiffCombinedResults.R)


### generating query signature from differentially expressed genes and KEGG pathways ppar and glyco
querySigDiffExp()
querySigKEGG()



###### run pharmacogenomics analysis
load("CMAP_gene_signatures.RData")

## 1- for the "diffExpAnalysis" query signature
library(PharmacoGx)
library(parallel)

ncor <- detectCores()
if (ncor>16) ncor <- 16

cl <- makeCluster(ncor)
load(paste("../Output/q850-diffExpAnalysis.RData", sep=""))
res <- parApply(CMAP.genePerturbations[,,c("tstat", "fdr")], 2, function(x, QSIG){
  source('./connectivityScore2.R')
  source('./combineTest.R')
  return(connectivityScore2(x=x,
                            y=QSIG,
                            method="gsea", nperm=1000))
}, cl = cl, QSIG=qq)

stopCluster(cl)
rownames(res) <- c("Connectivity", "P Value")
res <- t(res)
res <- res[order(res[,1], decreasing=T),]
save(list = ls(all.names = TRUE), file = paste("../Output/qSig",length(qq),"-diffExpAnalysis.RData", sep=""))

## 2- for the KEGG-based query signature
ncor <- detectCores()
if (ncor>16) ncor <- 16

cl <- makeCluster(ncor)

load(paste("../Output/q-KEGG.RData", sep=""))
res <- parApply(CMAP.genePerturbations[,,c("tstat", "fdr")], 2, function(x, QSIG){
  source('./connectivityScore2.R')
  source('./combineTest.R')
  return(connectivityScore2(x=x,
                            y=QSIG,
                            method="gsea", nperm=1000))
}, cl = cl, QSIG=qq)

stopCluster(cl)
rownames(res) <- c("Connectivity", "P Value")
res <- t(res)
res <- res[order(res[,1], decreasing=T),]
save(list = ls(all.names = TRUE), file = paste("../Output/qSig-KEGG-diffExpAnalysis.RData", sep=""))


## plot results
plotDiffCombinedResults()


