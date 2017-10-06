



rm(list=ls())
library(xlsx)

dat <- read.xlsx("../Data/BellySole.xlsx",1)


###NOTE:
## added on aug 14 ...
rownames(dat) = dat$NA.
dat = dat[,-1]



### transformation
dat <- log2(dat)
dat <- as.matrix(v)


x <- scale(dat)
x <- t(x)
library(gplots)

mydist=function(c) {as.dist(1-cor(t(c)))}
myclust=function(c) {hclust(c,method="complete")}


pdf("../Output/Fig1d-clustering.pdf")

color.map <- function(x) { if (strsplit(x, " ")[[1]][1] == "Belly") "#c51b8a"
  else if (strsplit(x, " ")[[1]][1] == "Sole") "orange"
  }


classColors <- unlist(lapply(colnames(x), color.map))

#
arrg <- heatmap.2(x, hclustfun=myclust, distfun = mydist, col=bluered, trace="none", ColSideColors = classColors, density.info = "none", Rowv= FALSE, labCol=FALSE, cexRow = 1, dendrogram="column") # , Rowv = FALSE, Colv = FALSE, dendogram = "none")
par(lend = 1)       
legend(x=0,y=0.75,    # location of the legend on the heatmap plot
       legend = c("Thin ECM","Dense ECM") ,
       
       col = c("#c51b8a","orange"),  # color key
       lty= 1,             # line style
       lwd = 10 ,           # line width
       cex = .7
)


dev.off()

