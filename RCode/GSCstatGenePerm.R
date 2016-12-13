## function adapted from Piano package: https://bioconductor.org/packages/release/bioc/html/piano.html


GSCstatGenePerm <- function (dummy, statistics, signs, gsc, statType, method, nGenes,
          nGenesUp, nGenesDn, nPerm, gseaParam) 
{
  source("./calcGeneSetStat.R")
  
  res <- list()
  res$gsStatsAllPerm <- list()
  res$gsStatsAllTestUpPerm <- list()
  res$gsStatsAllTestDnPerm <- list()
  res$gsStatsAbsPerm <- list()
  res$gsStatsUpPerm <- list()
  res$gsStatsDnPerm <- list()
  for (iContrast in 1:ncol(statistics)) {
    statsContrast <- statistics[, iContrast]
    if (statType %in% c("p-signed", "F-signed")) {
      signsContrast <- signs[, iContrast]
      statsContrastUp <- statsContrast[signsContrast > 
                                         0]
      statsContrastDn <- abs(statsContrast[signsContrast < 
                                             0])
    }
    else {
      statsContrastUp <- statsContrast[statsContrast > 
                                         0]
      statsContrastDn <- abs(statsContrast[statsContrast < 
                                             0])
    }
    if (method == "gsea") {
      statsContrastSorted <- sort(statistics[, iContrast], 
                                  decreasing = TRUE)
    }
    if (method %in% c("stouffer", "reporter", "tailStrength", 
                      "wilcoxon", "mean", "median", "sum") & statType == 
        "p-signed") {
      statsContrastTestUp <- abs(statsContrast)
      statsContrastTestUp[signsContrast > 0] <- statsContrastTestUp[signsContrast > 
                                                                      0]/2
      statsContrastTestUp[signsContrast < 0] <- 1 - statsContrastTestUp[signsContrast < 
                                                                          0]/2
      statsContrastTestUp[signsContrast == 0] <- NA
      statsContrastTestUp <- statsContrastTestUp[!is.na(statsContrastTestUp)]
      statsContrastTestDn <- 1 - statsContrastTestUp
      statsContrastTestUp[statsContrastTestUp < 1e-100] <- 1e-100
      statsContrastTestUp[statsContrastTestUp > 1] <- 1
      statsContrastTestDn[statsContrastTestDn < 1e-100] <- 1e-100
      statsContrastTestDn[statsContrastTestDn > 1] <- 1
    }
    gsSizes <- unique(nGenes[, iContrast])
    gsSizesUp <- unique(nGenesUp[, iContrast])
    gsSizesDn <- unique(nGenesDn[, iContrast])
    gsStatsAllMatrix <- matrix(nrow = length(gsSizes), ncol = nPerm)
    gsStatsAllTestUpMatrix <- matrix(nrow = length(gsSizes), 
                                     ncol = nPerm)
    gsStatsAllTestDnMatrix <- matrix(nrow = length(gsSizes), 
                                     ncol = nPerm)
    gsStatsAbsMatrix <- matrix(nrow = length(gsSizes), ncol = nPerm)
    gsStatsUpMatrix <- matrix(nrow = length(gsSizesUp), ncol = nPerm)
    gsStatsDnMatrix <- matrix(nrow = length(gsSizesDn), ncol = nPerm)
    for (iGeneSet in 1:length(gsSizes)) {
      gsStatsAll <- rep(NA, nPerm)
      gsStatsAllTestUp <- rep(NA, nPerm)
      gsStatsAllTestDn <- rep(NA, nPerm)
      gsStatsAbs <- rep(NA, nPerm)
      nGenesInSet <- gsSizes[iGeneSet]
      for (iPerm in 1:nPerm) {
        if (method == "fisher") {
          rInd <- sample(1:length(statsContrast), nGenesInSet)
          randStats <- statsContrast[rInd]
          gsStatsAbs[iPerm] <- calcGeneSetStat(abs(randStats), 
                                               method)
        }
        if (method %in% c("stouffer", "reporter", "tailStrength")) {
          rInd <- sample(1:length(statsContrast), nGenesInSet)
          randStats <- statsContrast[rInd]
          if (statType == "p-signed") {
            rInd <- sample(1:length(statsContrastTestUp), 
                           nGenesInSet)
            randStatsTestUp <- statsContrastTestUp[rInd]
            randStatsTestDn <- statsContrastTestDn[rInd]
            gsStatsAllTestUp[iPerm] <- calcGeneSetStat(randStatsTestUp, 
                                                       method)
            gsStatsAllTestDn[iPerm] <- calcGeneSetStat(randStatsTestDn, 
                                                       method)
          }
          gsStatsAbs[iPerm] <- calcGeneSetStat(abs(randStats), 
                                               method)
        }
        if (method == "wilcoxon") {
          rInd <- sample(1:length(statsContrast), nGenesInSet)
          randStats <- statsContrast[rInd]
          statsNotInSet <- statsContrast[!c(1:length(statsContrast)) %in% 
                                           rInd]
          if (statType == "p-signed") {
            rInd <- sample(1:length(statsContrastTestUp), 
                           nGenesInSet)
            randStatsTestUp <- statsContrastTestUp[rInd]
            randStatsTestDn <- statsContrastTestDn[rInd]
            statsNotInSetTestUp <- statsContrastTestUp[!c(1:length(statsContrastTestUp)) %in% 
                                                         rInd]
            statsNotInSetTestDn <- statsContrastTestDn[!c(1:length(statsContrastTestDn)) %in% 
                                                         rInd]
            gsStatsAllTestUp[iPerm] <- calcGeneSetStat(randStatsTestUp, 
                                                       "wilcoxon_fast", statsNotInSetTestUp)[1]
            gsStatsAllTestDn[iPerm] <- calcGeneSetStat(randStatsTestDn, 
                                                       "wilcoxon_fast", statsNotInSetTestDn)[1]
          }
          else if (statType == "t") {
            gsStatsAll[iPerm] <- calcGeneSetStat(randStats, 
                                                 "wilcoxon_fast", statsNotInSet)[1]
          }
          gsStatsAbs[iPerm] <- calcGeneSetStat(abs(randStats), 
                                               "wilcoxon_fast", abs(statsNotInSet))[1]
        }
        if (method %in% c("mean", "median", "sum")) {
          rInd <- sample(1:length(statsContrast), nGenesInSet)
          randStats <- statsContrast[rInd]
          if (statType == "p-signed") {
            rInd <- sample(1:length(statsContrastTestUp), 
                           nGenesInSet)
            randStatsTestUp <- statsContrastTestUp[rInd]
            randStatsTestDn <- statsContrastTestDn[rInd]
            gsStatsAllTestUp[iPerm] <- calcGeneSetStat(randStatsTestUp, 
                                                       method)
            gsStatsAllTestDn[iPerm] <- calcGeneSetStat(randStatsTestDn, 
                                                       method)
          }
          else if (statType == "t") {
            gsStatsAll[iPerm] <- calcGeneSetStat(randStats, 
                                                 method)
          }
          gsStatsAbs[iPerm] <- calcGeneSetStat(abs(randStats), 
                                               method)
        }
        if (method == "maxmean") {
          rInd <- sample(1:length(statsContrast), nGenesInSet)
          randStats <- statsContrast[rInd]
          gsStatsAbs[iPerm] <- calcGeneSetStat(randStats, 
                                               method)
        }
        if (method == "gsea") {
          rInd <- sample(1:length(statsContrast), nGenesInSet)
          gsStatsAll[iPerm] <- calcGeneSetStat(rInd, 
                                               method, statsContrastSorted, gseaParam)
        }
        if (method == "page") {
          rInd <- sample(1:length(statsContrast), nGenesInSet)
          randStats <- statsContrast[rInd]
          gsStatsAll[iPerm] <- calcGeneSetStat(randStats, 
                                               method, statsContrast)
        }
      }
      gsStatsAllMatrix[iGeneSet, ] <- gsStatsAll
      gsStatsAllTestUpMatrix[iGeneSet, ] <- gsStatsAllTestUp
      gsStatsAllTestDnMatrix[iGeneSet, ] <- gsStatsAllTestDn
      gsStatsAbsMatrix[iGeneSet, ] <- gsStatsAbs
    }
    if (statType %in% c("t", "p-signed", "F-signed") & method %in% 
        c("fisher", "stouffer", "reporter", "tailStrength", 
          "mean", "median", "sum", "wilcoxon")) {
      for (iGeneSet in 1:length(gsSizesUp)) {
        gsStatsUp <- rep(NA, nPerm)
        nGenesInSet <- gsSizesUp[iGeneSet]
        if (nGenesInSet > 0) {
          for (iPerm in 1:nPerm) {
            if (method %in% c("fisher", "stouffer", "reporter", 
                              "tailStrength", "mean", "median", "sum")) {
              rInd <- sample(1:length(statsContrastUp), 
                             nGenesInSet)
              randStats <- statsContrastUp[rInd]
              gsStatsUp[iPerm] <- calcGeneSetStat(randStats, 
                                                  method)
            }
            if (method == "wilcoxon") {
              rInd <- sample(1:length(statsContrastUp), 
                             nGenesInSet)
              randStats <- statsContrastUp[rInd]
              statsNotInSet <- statsContrastUp[!c(1:length(statsContrastUp)) %in% 
                                                 rInd]
              gsStatsUp[iPerm] <- calcGeneSetStat(randStats, 
                                                  "wilcoxon_fast", statsNotInSet)[1]
            }
          }
        }
        gsStatsUpMatrix[iGeneSet, ] <- gsStatsUp
      }
    }
    if (statType %in% c("t", "p-signed", "F-signed") & method %in% 
        c("fisher", "stouffer", "reporter", "tailStrength", 
          "mean", "median", "sum", "wilcoxon")) {
      for (iGeneSet in 1:length(gsSizesDn)) {
        gsStatsDn <- rep(NA, nPerm)
        nGenesInSet <- gsSizesDn[iGeneSet]
        if (nGenesInSet > 0) {
          for (iPerm in 1:nPerm) {
            if (method %in% c("fisher", "stouffer", "reporter", 
                              "tailStrength", "mean", "median", "sum")) {
              rInd <- sample(1:length(statsContrastDn), 
                             nGenesInSet)
              randStats <- statsContrastDn[rInd]
              gsStatsDn[iPerm] <- calcGeneSetStat(randStats, 
                                                  method)
            }
            if (method == "wilcoxon") {
              rInd <- sample(1:length(statsContrastDn), 
                             nGenesInSet)
              randStats <- statsContrastDn[rInd]
              statsNotInSet <- statsContrastDn[!c(1:length(statsContrastDn)) %in% 
                                                 rInd]
              gsStatsDn[iPerm] <- calcGeneSetStat(randStats, 
                                                  "wilcoxon_fast", statsNotInSet)[1]
            }
          }
        }
        gsStatsDnMatrix[iGeneSet, ] <- gsStatsDn
      }
    }
    rownames(gsStatsAllMatrix) <- gsSizes
    rownames(gsStatsAllTestUpMatrix) <- gsSizes
    rownames(gsStatsAllTestDnMatrix) <- gsSizes
    rownames(gsStatsAbsMatrix) <- gsSizes
    rownames(gsStatsUpMatrix) <- gsSizesUp
    rownames(gsStatsDnMatrix) <- gsSizesDn
    res$gsStatsAllPerm[[iContrast]] <- gsStatsAllMatrix
    res$gsStatsAllTestUpPerm[[iContrast]] <- gsStatsAllTestUpMatrix
    res$gsStatsAllTestDnPerm[[iContrast]] <- gsStatsAllTestDnMatrix
    res$gsStatsAbsPerm[[iContrast]] <- gsStatsAbsMatrix
    res$gsStatsUpPerm[[iContrast]] <- gsStatsUpMatrix
    res$gsStatsDnPerm[[iContrast]] <- gsStatsDnMatrix
  }
  return(res)
}