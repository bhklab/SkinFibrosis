

## this code clusters the temporal expression data (only fibrostic/irradiated group)
## then performs overrepresentation analysis on the two largest clusters ...


rm(list=ls())

library(piano)
library(xlsx)
library(WriteXLS)
library(GSA)
library(maSigPro)


## STEP 1: CLUSTERING
edesign <- matrix(nrow=9,ncol=3)
rownames(edesign) <- c("Radiated_0d_1", "Radiated_0d_2", "Radiated_0d_3" ,
                       "Radiated_6d_1", "Radiated_6d_2", "Radiated_6d_3", 
                       "Radiated_20d_1", "Radiated_20d_2", "Radiated_20d_3"
)


colnames(edesign) <- c("Time", "Replicates", "Radiated")
edesign[,"Time"] <- rep( c(0,6,20),each=3)
edesign[,"Replicates"] <- rep( c(1:3),each=3) 
edesign[,"Radiated"] <- c(rep(1,9))

exprs <- as.matrix(read.table("../Data/genes.fpkm_table", header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))

## re-order according to the design matrix rownames
exprs3 <- matrix(nrow=dim(exprs)[1], ncol=9)
exprs3[,1] <- exprs[,15] ### N  --> N_0
exprs3[,2] <- exprs[,13] # N_1
exprs3[,3] <- exprs[,14] # N_2
#
exprs3[,4] <- exprs[,8] ## S --> S_0
exprs3[,5] <- exprs[,9] #S_1
exprs3[,6] <- exprs[,7] #S_2
#
exprs3[,7:9] <- exprs[,1:3] ## L_0, L_1, L_2
# 
colnames(exprs3) <- rownames(edesign) 
rownames(exprs3) <- rownames(exprs)

exprs3 <- log2(exprs3+1)
exprs3 <- exprs3[apply(exprs3[,-1], 1, function(x) !all(x==0)),]


## only include genes in the union of diffrerentially expressed genes 
d <- read.table("../Data/gene_exp.diff")
d2 <- subset(d, d$V5=='L' & d$V6 =='S' & d$V11 !="-nan") ## NOTE
d22 <- subset(d, d$V5=='N' & d$V6 =='S' & d$V11 !="-nan") ## NOTE

## remove KRT
indices <- gdata::startsWith(as.character(d2$V1), "Krt")
sum(indices)
#[1] 107
d2 <- d2[ !gdata::startsWith(as.character(d2$V1), "Krt"),]

indices <- gdata::startsWith(as.character(d22$V1), "Krt")
sum(indices)
#[1] 107
d22 <- d22[ !gdata::startsWith(as.character(d22$V1), "Krt"),]


d2$V13 <- as.numeric(as.character(d2$V13))
d2$V10 <- as.numeric(as.character(d2$V10))
d22$V13 <- as.numeric(as.character(d22$V13))
d22$V10 <- as.numeric(as.character(d22$V10))

## subset the most significant genes
d2 <- d2[d2$V13<0.05,]
d22 <- d22[d22$V13<0.05,]

dAllSig <- rbind(d2, d22)
dim(dAllSig)
##[1] 3898   14

length(unique(dAllSig$V1))
## 3500
genes <- unique(dAllSig$V1) # --> 3500 genes  

######### filter out the exprs3 
exprs3 <- exprs3[rownames(exprs3) %in% genes,]
dim(exprs3)
# 3500 
#########

## apply maSigPro pipeline
pdf("../Output/Fig1/clusteringResults-KRTFiltered-fibrosisOnly.pdf")
x <- see.genes(exprs3, edesign=edesign, newX11 =  FALSE)
dev.off()  
## 
cluster1 <- exprs3[x$cut==1,]
cluster2 <- exprs3[x$cut==2,]
cluster3 <- exprs3[x$cut==3,]
cluster4 <- exprs3[x$cut==4,] 
cluster5 <- exprs3[x$cut==5,]
cluster6 <- exprs3[x$cut==6,]
cluster7 <- exprs3[x$cut==7,]
cluster8 <- exprs3[x$cut==8,]
cluster9 <- exprs3[x$cut==9,]
save(list = ls(all.names = TRUE), file = "../Output/Fig1/clusteringResults-KRTFiltered-fibrosisOnly.RData")


## STEP 2: OVER-REP. ANALYSIS
## pathways have been downloaded from consensusPAthDB website
paths <- read.xlsx("../Data/MouseReactome.xlsx",1)

## subset to kegg only
#paths <- paths[paths$source=="KEGG",]
#paths <- paths[paths$source=="MouseCyc",]

paths$pathway <- as.character(paths$pathway)
paths$hgnc_symbol_ids <- toupper(as.character(paths$hgnc_symbol_ids))
gTogs <- data.frame()
for (i in 1:dim(paths)[1]) {
  x <- unlist(strsplit(paths$hgnc_symbol_ids[i], ","))
  gTogs <- rbind(gTogs,(cbind(x, GSetName=paste(paths$pathway[i], paths$source[i], sep="%"))))
}
length(unique(gTogs$GSetName))
a <- loadGSC(gTogs)

d <- read.table("../Data/gene_exp.diff")
d2 <- subset(d, d$V5=='L' & d$V6 =='S' & d$V11 !="-nan") ## NOTE
d22 <- subset(d, d$V5=='N' & d$V6 =='S' & d$V11 !="-nan") ## NOTE


### OVER REP for two largest clusters ...
## cluster 4
universe <- toupper(intersect(d2$V1, d22$V1))
genes <- toupper(rownames(cluster4))
significant <- rep(1, length(universe))
significant[which(universe %in% genes)] <- 0
res <- runGSAhyper(genes=universe, significant,  gsc=a, adjMethod="fdr")
result <- res$resTab[order(res$resTab[, "p-value"]),]
Sys.setlocale('LC_ALL','C') 
write.xlsx2(result, file= paste("../Output/Fig1/cluster4-ReactomeOnly-fig1h.xlsx", sep=""))
## cluster 2
universe <- toupper(intersect(d2$V1, d22$V1))
genes <- toupper(rownames(cluster2))
significant <- rep(1, length(universe))
significant[which(universe %in% genes)] <- 0
res <- runGSAhyper(genes=universe, significant,  gsc=a, adjMethod="fdr")
result <- res$resTab[order(res$resTab[, "p-value"]),]
Sys.setlocale('LC_ALL','C') 
write.xlsx2(result, file= paste("../Output/Fig1/cluster2-ReactomeOnly-fig1h.xlsx", sep=""))





