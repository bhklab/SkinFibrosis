

querySigDiffExp <- function() {

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


prbNames <- colnames(ndata.all) ## after removing the duplicated probe ids ...
ndata.all <- read.csv("../Data/modTTESTcorALLpvalues_alldatafiltered.txt", sep="\t", comment.char="#", header=TRUE) # 33788    80
rownames(ndata.all) <- paste("probe", ndata.all[ , "ProbeID"], sep=".")
ndata.all <- ndata.all[ prbNames ,c("ProbeID", "p..Corr.", "p", "Regulation", "FC..abs.", "FC", "Log.FC", "Symbol", "Entrez_Gene_ID")]

tt <- ndata.all
tt <- tt[tt$Entrez_Gene_ID!="",]
tt <- tt[!is.na(tt$Entrez_Gene_ID),]

## remove low fdr
tt <- tt[(tt[, "p..Corr."]<0.05),]

###  SORT tt according to FC
tt <- tt[order(tt$FC, decreasing = TRUE),]

tt2 <- tt
######
q <- tt2$FC
names(q) <- rownames(tt2)
xxentr <- lapply(tt2$Entrez_Gene_ID,function(x) paste("geneid.",x, sep=''))
names(q) <- xxentr

sz <- 425
qup <- head(q, sz)
# reversing the sign
qup <- unlist(lapply(qup, function(x) {if (x>0) {x=-1} else {x= 1}} ))
qdn <- tail(q, sz)
# reversing the sign
qdn <- unlist(lapply(qdn, function(x) {if (x>0) {x=-1} else {x= 1}} ))
qq <- c(qup,qdn)
save(qq, file=paste("./Output/q", (sz)*2, "-diffExpAnalysis.RData", sep=""))


}





