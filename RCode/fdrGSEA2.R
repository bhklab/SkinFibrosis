## function adapted from Piano package: https://bioconductor.org/packages/release/bioc/html/piano.html


fdrGSEA2 <- function (gsStatsAll, gsStatsAllPerm, nGenes, signMethod)
{
  res <- list()
  pValuesAllUpAdj <- vector()
  pValuesAllDnAdj <- vector()
  for (iContrast in 1:ncol(gsStatsAll)) {
    randBgMat <- gsStatsAllPerm[[iContrast]]
    randBgMatNorm <- randBgMat
    for (iGeneSetSize in 1:nrow(randBgMat)) {
      geneSetBg <- randBgMat[iGeneSetSize, ]
      posMean <- max(c(mean(geneSetBg[geneSetBg >= 0]), 
                       0), na.rm = TRUE)
      negMean <- max(c(abs(mean(geneSetBg[geneSetBg <= 
                                            0])), 0), na.rm = TRUE)
      geneSetBg[geneSetBg > 0] <- geneSetBg[geneSetBg > 
                                              0]/posMean
      geneSetBg[geneSetBg < 0] <- geneSetBg[geneSetBg < 
                                              0]/negMean
      randBgMatNorm[iGeneSetSize, ] <- geneSetBg
    }
    NES <- rep(NA, nrow(gsStatsAll))
    randBgMatNormFull <- as.data.frame(matrix(nrow = nrow(gsStatsAll), 
                                              ncol = ncol(randBgMat)))
    for (iGeneSet in 1:nrow(gsStatsAll)) {
      if (signMethod == "geneperm") {
        sizeGS <- nGenes[iGeneSet, iContrast]
        randBgMatNormFull[iGeneSet, ] <- randBgMatNorm[as.character(sizeGS), 
                                                       ]
      }
      else {
        randBgMatNormFull[iGeneSet, ] <- randBgMatNorm[iGeneSet, 
                                                       ]
      }
      ES <- gsStatsAll[iGeneSet, iContrast]
      if (signMethod == "geneperm") {
        tmp <- randBgMat[as.character(sizeGS), ]
      }
      else {
        tmp <- randBgMat[iGeneSet, ]
      }
      if (ES > 0) {
        NES[iGeneSet] <- ES/max(c(mean(tmp[tmp >= 0]), 
                                  0), na.rm = TRUE)
      }
      else if (ES < 0) {
        NES[iGeneSet] <- ES/max(c(abs(mean(tmp[tmp <= 
                                                 0])), 0), na.rm = TRUE)
      }
      else {
        NES[iGeneSet] <- 0
      }
    }
    FDRup <- rep(NA, nrow(gsStatsAll))
    FDRdn <- rep(NA, nrow(gsStatsAll))
    for (iGeneSet in 1:nrow(gsStatsAll)) {
      if (NES[iGeneSet] > 0) {
        tmp1 <- sum(randBgMatNormFull >= NES[iGeneSet])/sum(randBgMatNormFull >= 
                                                              0)
        tmp2 <- sum(NES >= NES[iGeneSet])/sum(NES >= 
                                                0)
        FDRup[iGeneSet] <- tmp1/tmp2
        FDRdn[iGeneSet] <- NA
      }
      else if (NES[iGeneSet] < 0) {
        tmp1 <- sum(randBgMatNormFull <= NES[iGeneSet])/sum(randBgMatNormFull <= 
                                                              0)
        tmp2 <- sum(NES <= NES[iGeneSet])/sum(NES <= 
                                                0)
        FDRdn[iGeneSet] <- tmp1/tmp2
        FDRup[iGeneSet] <- NA
      }
      else {
        FDRup[iGeneSet] <- 1
        FDRdn[iGeneSet] <- 1
      }
      if (max(c(FDRup[iGeneSet], 1), na.rm = TRUE) > 1) 
        FDRup[iGeneSet] <- 1
      if (max(c(FDRdn[iGeneSet], 1), na.rm = TRUE) > 1) 
        FDRdn[iGeneSet] <- 1
    }
    pValuesAllUpAdj <- cbind(pValuesAllUpAdj, FDRup)
    pValuesAllDnAdj <- cbind(pValuesAllDnAdj, FDRdn)
  }
  res$pValuesAllUpAdj <- pValuesAllUpAdj
  res$pValuesAllDnAdj <- pValuesAllDnAdj
  res$NES <- NES
  return(res)
}
