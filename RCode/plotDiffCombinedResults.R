


plotDiffCombinedResults <- function() { 

load("../Output/qSig850-diffExpAnalysis.RData")
iix <- which(rownames(res)=="caffeic acid")
resdiff <- res

load("../Output/qSig-KEGG.RData")
which(rownames(res)=="caffeic acid")
# # #[1] 307
resCombined <- res
# 

ii <- match(rownames(resdiff), rownames(resCombined))
all(rownames(resCombined)[ii] == rownames(resdiff)) 
resCombined <- resCombined[ii,]

comb <- cbind(resdiff[,"Connectivity"], resCombined[,"Connectivity"])
colnames(comb) <- c("resDiff","resCombined")

iix <- which(rownames(comb)=="caffeic acid")

sum(comb[,1]>0) ## 267
sum(comb[,2]>0) ## 491    
i1 <- comb[,1]>0
i2 <- comb[,2]>0
ix <- i1 & i2 

topdiffandcombined <- comb[ix,]

######
rownames(comb)[iix]
cafacid <- comb[iix,]
comb <- comb[-iix,]
pdf("../Output/resDiffCombined-qSig850-ppar-glyco.pdf")
plot(comb, xlab = "conn score (human diff exps)", ylab = "conn score (ppar-glyco)", ylim=c(-0.6,0.6), xlim=c(-0.6,0.6))

dev.off()
}





