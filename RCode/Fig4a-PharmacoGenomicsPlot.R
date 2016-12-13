

## This script generates the heatmap in figure 4a. from 6 ranked lists of drugs. The lists are results of running 
## the computational pharmacogenomics pipeline impelemented in : "PharmacogenomicsAnalysis.R"


rm(list=ls())

library(xlsx)

f <- list.files("../Data/drugResultsLists")
lst <- read.xlsx(file=f[1],1)  
rankMat <- matrix(nrow=dim(lst)[1], ncol=6)
rankMat[,1] <- order(1-lst$Connectivity)
rownames(rankMat) = lst$DrugName

x2 <- read.xlsx(file=f[2],1)  
x3 <- read.xlsx(file=f[3],1)  
x4 <- read.xlsx(file=f[4],1)  
x5 <- read.xlsx(file=f[5],1)  
x6 <- read.xlsx(file=f[6],1)  

for (i in rownames(rankMat)) {
  ## find the indices
  rankMat[i,2] <- which(i == x2$DrugName)
  rankMat[i,3] <- which(i == x3$DrugName)
  rankMat[i,4] <- which(i == x4$DrugName)
  rankMat[i,5] <- which(i == x5$DrugName)
  rankMat[i,6] <- which(i == x6$DrugName)
}

## each drug has 6 coordinates
rankSum <- rowSums(rankMat)
names(rankSum) <- rownames(rankMat)
## sort the rank sum  
xx <- rankMat[order(rankSum),]
xx[xx>100] <- NA

colnames(xx) <- c("Query Size = 150","Query Size = 200","Query Size = 250","Query Size = 300","Query Size = 350","Query Size = 400")

library(pheatmap)

## capitalize the drug names
simpleCap <- function(x) {
  substring(x, 1,1) <- toupper(substring(x, 1,1))
  return(x)
}
rownames(xx) <- lapply(rownames(xx), simpleCap)


pdf(paste("Fig4a-pharmaCoPlot.pdf",sep=""), width=8,height=10)
par(mar=c(10,16,4,2)+0.1,mgp=c(13,1,0))
pheatmap(xx[1:20,], cluster_rows = F,cluster_cols = F,legend = T, cellheight = 15,border_color = "black")
dev.off()


