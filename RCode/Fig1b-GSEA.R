

rm(list=ls())
library(snow)
library(parallel)
library(GSA)
source("runGSA2.R")
source("fdrGSEA2.R")
source("checkLoadArg.R")
source("GSCstatBatch.R")
source("GSCsignificanceBatch.R")
source("GSCstatGenePerm.R")
source("pvalFromFractionGenePerm.R")
source("loadGSC.R")


## list of characters to be removed from row and column names
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"


## read combined pheno data
pdata.all <- read.csv("../Data/SampleAnnotation.txt", sep="\t")
rownames(pdata.all) <- gsub(badchars, "_", pdata.all[ , "Sample.ID"])
## read combined data
ndata.all <- read.csv("../Data/modTTESTcorALLpvalues_alldatafiltered.txt", sep="\t", comment.char="#", header=TRUE) # 33788    80
annot <- ndata.all[ , c("Symbol", "Entrez_Gene_ID", "Definition", "Ontology_Process", "Ontology_Component", "Ontology_Function", "Chromosome", "Chromosome.Start.Index.Avadis.", "Chromosome.End.Index.Avadis.", "Probe_Sequence")]
rownames(annot) <- rownames(ndata.all) <- paste("probe", ndata.all[ , "ProbeID"], sep=".")

## select only normalized data
ndata.all <- ndata.all[ , grep(".normalized.", colnames(ndata.all)), drop=FALSE]
cc <- gsub(badchars, "_", sapply(strsplit(gsub("\\[", "", colnames(ndata.all)), split=","), function (x) { return (x[1]) }))
cc <- unlist(lapply(cc, function(x) gsub("X_", "", x)))
cc <- unlist(lapply(cc, function(x) gsub("__R_", "", x)))
cc <- unlist(lapply(cc, function(x) gsub("__N_", "", x)))

colnames(ndata.all) <-  unlist(lapply(cc, function(x) gsub("_normalized_", "", x)))
ndata.all <- t(ndata.all)
nn <- dimnames(ndata.all)
ndata.all <- apply(ndata.all, 2, as.numeric)
dimnames(ndata.all) <- nn
## reorder patients
sn <- intersect(rownames(pdata.all), rownames(ndata.all))
ndata.all <- ndata.all[sn, , drop=FALSE]
pdata.all <- pdata.all[sn, , drop=FALSE]


ndata.all <- ndata.all[!rownames(ndata.all) %in% c("P9_N","P15_N"), , drop=FALSE]
pdata.all <- pdata.all[!rownames(pdata.all) %in% c("P9_N","P15_N"), , drop=FALSE]


### NOTE!!!! added on Nov 14: filter out duplicated probes...
all(rownames(annot) == colnames(ndata.all))
annotSym <- annot$Symbol[annot$Symbol!=""]
ix <- duplicated(annotSym)
sum(ix)
dupSyms <- unique(annotSym[ix])

dels <- vector(mode="character")
for (i in 1:length(dupSyms)) {
  ii <- which(annot$Symbol %in% dupSyms[i])
  data <- ndata.all[,ii]
  x <- apply(data, 2, function(x) {IQR(x)})
  dels <- c(dels, names(x[x!=max(x)])) 
  #ndata.all <- ndata.all[, -c(which(colnames(ndata.all) %in% dels))]
}
## 
length(unique(dels)) == length(dels)
ndata.all <- ndata.all[, -c(which(colnames(ndata.all) %in% dels))]
dim(ndata.all)
#[1]    22 26383


## radiated 13
ndata.r <- ndata.all[grep("_R", rownames(ndata.all)), colnames(ndata.all), drop=FALSE]
pdata.r <- pdata.all[grep("_R", rownames(pdata.all)), , drop=FALSE]

## normal 11
ndata.n <- ndata.all[grep("_N", rownames(ndata.all)), colnames(ndata.all), drop=FALSE]
pdata.n <- pdata.all[grep("_N", rownames(pdata.all)), , drop=FALSE]


## diff expression
ndata.n <- t(ndata.n)
ndata.r <- t(ndata.r)


ndata.r <- 2 ^ ndata.r
ndata.n <- 2 ^ ndata.n
tt <- matrix(nrow=dim(ndata.n)[1], ncol=2 )
rownames(tt) <- rownames(ndata.n)
colnames(tt) <- c("logFC", "FC")

for(g in 1:dim(tt)[1]) {
  tt[g,"logFC"] <- log2(mean(ndata.r[g,])) - log2(mean(ndata.n[g,]))
  tt[g,"FC"] <- mean(ndata.r[g,]) / mean(ndata.n[g,])
}

#all(rownames(annot)==rownames(tt))
tt <- as.data.frame(tt)
tt$Symbol <- annot[rownames(annot) %in% rownames(tt), "Symbol"]
tt <- tt[tt$Symbol!="",]
tt <- tt[!is.na(tt$Symbol),]
dim(tt)


###### ###### ###### ###### ###### GSEA for GO TERMS ###### ###### ###### ###### ###### ###### 
for (sourceName in c("bp","mf","cc")) {
 gSets <- GSA.read.gmt(paste("c5.", sourceName, ".v5.2.symbols.gmt", sep=""))

 dfgSets <- as.data.frame(cbind(unlist(gSets$genesets))) ## genes
 dfgSNames <- as.data.frame(cbind(unlist(gSets$geneset.names))) 

 listS <- lapply(gSets$genesets,length)

 ##  now need to combine the df$V1 and dd
 gTogs <- data.frame(dfgSets$V1,rep(dfgSNames$V1,listS))
 ## filter out the empty gene names
 names(gTogs) <- c("V1","V2")
 gTogs <- gTogs[gTogs$V1!="",]

 gTogs$V1 <- toupper(gTogs$V1)
 gTogs$V1 <- gsub(badchars, "", gTogs$V1)
 a <- loadGSC(gTogs)
 ################################
 gg <- tt$logFC
 names(gg) <- tt$Symbol
 names(gg) <- gsub(badchars, "", names(gg))
 names(gg) <- toupper(names(gg))

 res <- runGSA2(geneLevelStats=gg, geneSetStat="gsea", gsc=a, nPerm = 1000, ncpus=10)

 gs <- piano::GSAsummaryTable(res) 
 gs$NES <- res$NES

 gs2 <- cbind(gs, "p"=NA, "p adj"=NA)
 gs2[ , "p"] <- gs2[ , "p (dist.dir.up)"]
 naix <- is.na(gs2[, "p"])
 gs2[naix, "p"] <- gs2[naix, "p (dist.dir.dn)"]
 gs2[ , "p adj"] <- gs2[ , "p adj (dist.dir.up)"]
 naix <- is.na(gs2[, "p adj"])
 gs2[naix, "p adj"] <- gs2[naix, "p adj (dist.dir.dn)"]
 gs2 <- gs2[!is.element(colnames(gs2), c("p (dist.dir.up)", "p (dist.dir.dn)", "p adj (dist.dir.up)", "p adj (dist.dir.dn)"))]
 gs2<- gs2[order(gs2[ , "p"]), ]

 zpix <- gs2$p==0
 nPerm=1000;
 gs2[zpix,"p"] = 1/(nPerm+1)  
 gs2 <- cbind(gs2, "fdr"=p.adjust(gs2[ , "p"], method="fdr"))

 ## save all workspace 
 save(res, gs2, gs, file = paste("../Output/GSEA-RvsN-Human-Fibrosis-", sourceName, ".RData", sep=""))
 write.xlsx(gs2, file = paste("../Output/GSEA-RvsN-Human-Fibrosis-", sourceName, ".xlsx", sep=""))
}



###### ###### ###### ###### ###### GSEA for humanCyc and Reactome and Kegg ###### ###### ###### ###### ###### ###### 
aa <- read.csv("../Data/CPDB_pathways_genes-human.tab", sep="\t", comment.char="#")
## subset to kegg only
for (sourcName in c("humanCyc", "Reactome", "KEGG")) {
 ## 328 HumanCyc
 ## 290 KEGG
 ##  1588   Reactome
 aa <- aa[aa$source==sourceName,] 
 aa$pathway <- as.character(aa$pathway)
 aa$hgnc_symbol_ids <- as.character(aa$hgnc_symbol_ids)
 aa$hgnc_symbol_ids <- toupper(aa$hgnc_symbol_ids)
 gTogs <- data.frame()
 for (i in 1:dim(aa)[1]) {
   x <- unlist(strsplit(aa$hgnc_symbol_ids[i], ","))
   x <- gsub(badchars, "", x)
   gTogs <- rbind(gTogs,(cbind(x, GSetName=paste(aa$pathway[i], aa$source[i], sep="%"))))
 }
 length(unique(gTogs$GSetName))
 a <- loadGSC(gTogs)

 gg <- tt$logFC
 names(gg) <- tt$Symbol
 names(gg) <- gsub(badchars, "", names(gg))
 names(gg) <- toupper(names(gg))

 res <- runGSA2(geneLevelStats=gg, geneSetStat="gsea", gsc=a, nPerm = 1000, ncpus = 10)

 gs <- piano::GSAsummaryTable(res)  #, save=TRUE), file="gsaResTab-GO-mf-Oct27.xls")
 gs$NES <- res$NES

 gs2 <- cbind(gs, "p"=NA, "p adj"=NA)
 gs2[ , "p"] <- gs2[ , "p (dist.dir.up)"]
 naix <- is.na(gs2[, "p"])
 gs2[naix, "p"] <- gs2[naix, "p (dist.dir.dn)"]
 gs2[ , "p adj"] <- gs2[ , "p adj (dist.dir.up)"]
 naix <- is.na(gs2[, "p adj"])
 gs2[naix, "p adj"] <- gs2[naix, "p adj (dist.dir.dn)"]
 gs2 <- gs2[!is.element(colnames(gs2), c("p (dist.dir.up)", "p (dist.dir.dn)", "p adj (dist.dir.up)", "p adj (dist.dir.dn)"))]
 gs2<- gs2[order(gs2[ , "p"]), ]

 zpix <- gs2$p==0
 nPerm=1000;
 gs2[zpix,"p"] = 1/(nPerm+1)  
 gs2 <- cbind(gs2, "fdr"=p.adjust(gs2[ , "p"], method="fdr"))

 ## save all workspace 
 save(res, gs2, gs, file = paste("../Output/GSEA-RvsN-Human-Fibrosis-", sourceName, ".RData", sep=""))
 write.xlsx(gs2, file = paste("../Output/GSEA-RvsN-Human-Fibrosis-", sourceName, ".xlsx", sep=""))
 
}


