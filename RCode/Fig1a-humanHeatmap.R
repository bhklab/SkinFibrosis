

## This script generates a heatmap from human data for fig1a


rm(list=ls())

## list of characters to be removed from row and column names
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"


## read annotation data
pdata.all <- read.csv("../Data/SampleAnnotation.txt", sep="\t")
rownames(pdata.all) <- gsub(badchars, "_", pdata.all[ , "Sample.ID"])

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

## remove the outliers
ndata.all <- ndata.all[!rownames(ndata.all) %in% c("P9_N","P15_N"), , drop=FALSE]
pdata.all <- pdata.all[!rownames(pdata.all) %in% c("P9_N","P15_N"), , drop=FALSE]


###  filter out duplicated probes...
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
}
## 
length(unique(dels)) == length(dels)
ndata.all <- ndata.all[, -c(which(colnames(ndata.all) %in% dels))]
dim(ndata.all)
#[1]    22 26383


rownames(ndata.all) <- unlist(lapply(rownames(ndata.all), function(x) {unlist(strsplit(x, "_"))[2]}))


## NOTE: pre-processing of the genes before performing clustering 
## find the most variant genes and pick the top 1000
vv <- apply(ndata.all, 2, function(x) {return(var(x))})
vv <- sort(vv, decreasing = TRUE)
top1000 <- vv[1:1000]

## mean center each gene before clustering ...
ww <- scale(ndata.all)
ww <- ww[,colnames(ww) %in% names(top1000)]

x <- t(ww)
library(gplots)


mydist=function(c) {as.dist(1-cor(t(c)))}
myclust=function(c) {hclust(c,method="complete")}


pdf("../Output/Fig1/Fig1a-HumanHeatmap.pdf")
par(mar=c(1,3,1,2)+0.1,mgp=c(3,1,0))
color.map <- function(x) { if (x=="R") "#c51b8a" else "#31a354" }
classColors <- unlist(lapply(colnames(x), color.map))
arrg <- heatmap.2(x, hclustfun=myclust, distfun=mydist, col=bluered, trace="none", ColSideColors = classColors, labRow=FALSE, density.info = "none", dendrogram = "column") # , Rowv = FALSE, Colv = FALSE, dendogram = "none")
par(lend = 1)       
legend(x=0.05,y=0.75,    # location of the legend on the heatmap plot
       legend = c("R", "N") ,
       col = c("#c51b8a","#31a354"),  # color key
       lty= 1,             # line style
       lwd = 10 ,           # line width
       cex = .7
)


dev.off()



