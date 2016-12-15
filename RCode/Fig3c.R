

rm(list=ls())

exprs <- as.matrix(read.table("../Data/genes.fpkm_table", header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
exprs <- log2(exprs+1)
exprs <- exprs[,c("L_0","L_1","L_2","K_0","K_1","K_2")]
exprs <- t(exprs)


## filter out the top 1000 most variant genes ...
## NOTE: pre-processing of the genes before doing clustering 
## find the most variant genes and pick the top 1000
#vv <- var(ndata.all)
vv <- apply(exprs, 2, function(x) {return(var(x))})
vv <- sort(vv, decreasing = TRUE)
top1000 <- vv[1:1000]

## mean center each gene before clustering ...
ww <- scale(exprs)
ww <- ww[,colnames(ww) %in% names(top1000)]


x <- t(ww)
library(gplots)


mydist=function(c) {as.dist(1-cor(t(c)))}
myclust=function(c) {hclust(c,method="complete")}

pdf("../Output/ADSCvsR-Heatmap.pdf")
#pdf("test.pdf")
par(mar=c(1,3,1,2)+0.1,mgp=c(3,1,0))
color.map <- function(x) { if (x %in% c("L_0", "L_1", "L_2")) "#c51b8a" else "#31a354" }
classColors <- unlist(lapply(colnames(x), color.map))
arrg <- heatmap.2(x, hclustfun=myclust, distfun=mydist, col=bluered, trace="none", ColSideColors = classColors, labRow=FALSE, density.info = "none", dendrogram = "column") # , Rowv = FALSE, Colv = FALSE, dendogram = "none")
par(lend = 1)       
legend(x=0.05,y=0.75,    # location of the legend on the heatmap plot
       legend = c("L", "K") ,
       col = c("#c51b8a", "#31a354"),  # color key
       lty= 1,             # line style
       lwd = 10,           # line width
       cex = .7
)

dev.off()


