
## function adapted from Piano package: https://bioconductor.org/packages/release/bioc/html/piano.html


GSCstatBatch <- function (statistics, statType, gsc, method, signMethod, gseaParam,
          signs) 
{
  source("calcGeneSetStat.R")
  
  res <- list()
  res$statsAll <- vector()
  res$statsAllTestUp <- vector()
  res$statsAllTestDn <- vector()
  res$statsAbs <- vector()
  res$statsUp <- vector()
  res$statsDn <- vector()
  res$nGenes <- vector()
  res$nGenesUp <- vector()
  res$nGenesDn <- vector()
  res$pValuesAll <- vector()
  res$pValuesAllUp <- vector()
  res$pValuesAllDn <- vector()
  res$pValuesAbs <- vector()
  res$pValuesUp <- vector()
  res$pValuesDn <- vector()
  for (iContrast in 1:ncol(statistics)) {
    gsStatsAll <- rep(NA, length(gsc))
    gsStatsAllTestUp <- rep(NA, length(gsc))
    gsStatsAllTestDn <- rep(NA, length(gsc))
    gsStatsAbs <- rep(NA, length(gsc))
    gsStatsUp <- rep(NA, length(gsc))
    gsStatsDn <- rep(NA, length(gsc))
    nGenes <- rep(0, length(gsc))
    nGenesUp <- rep(0, length(gsc))
    nGenesDn <- rep(0, length(gsc))
    pValuesAll <- rep(NA, length(gsc))
    pValuesAllUp <- rep(NA, length(gsc))
    pValuesAllDn <- rep(NA, length(gsc))
    pValuesAbs <- rep(NA, length(gsc))
    pValuesUp <- rep(NA, length(gsc))
    pValuesDn <- rep(NA, length(gsc))
    statsContrast <- statistics[, iContrast]
    if (statType %in% c("p-signed", "F-signed")) {
      signsContrast <- signs[, iContrast]
    }
    if (statType == "t") {
      signsContrast <- sign(statsContrast)
    }
    if (method %in% c("stouffer", "reporter", "tailStrength", 
                      "wilcoxon", "wilcoxon_fast", "mean", "median", "sum") & 
        statType == "p-signed") {
      statsContrastTestUp <- abs(statsContrast)
      statsContrastTestUp[signsContrast > 0] <- statsContrastTestUp[signsContrast > 
                                                                      0]/2
      statsContrastTestUp[signsContrast < 0] <- 1 - statsContrastTestUp[signsContrast < 
                                                                          0]/2
      statsContrastTestUp[signsContrast == 0] <- NA
      statsContrastTestDn <- 1 - statsContrastTestUp
      statsContrastTestUp[statsContrastTestUp < 1e-100] <- 1e-100
      statsContrastTestUp[statsContrastTestUp > 1] <- 1
      statsContrastTestDn[statsContrastTestDn < 1e-100] <- 1e-100
      statsContrastTestDn[statsContrastTestDn > 1] <- 1
    }
    for (iGeneSet in 1:length(gsc)) {
      if (method == "fisher") {
        indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
        statsGenesInSet <- statsContrast[indGenesInSet]
        nGenes[iGeneSet] <- length(statsGenesInSet)
        if (nGenes[iGeneSet] > 0) {
          if (exists("signsContrast")) {
            nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet] > 
                                        0)
            nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet] < 
                                        0)
          }
          gsStatsAbs[iGeneSet] <- calcGeneSetStat(abs(statsGenesInSet),
                                                  method)
          if (statType == "p-signed" & signMethod != 
              "sampleperm") {
            signsGenesInSet <- signsContrast[indGenesInSet]
            sel <- statsGenesInSet[signsGenesInSet > 
                                     0]
            nGenesUp[iGeneSet] <- length(sel)
            if (nGenesUp[iGeneSet] > 0) {
              gsStatsUp[iGeneSet] <- calcGeneSetStat(sel,
                                                     method)
            }
            sel <- abs(statsGenesInSet[signsGenesInSet < 
                                         0])
            nGenesDn[iGeneSet] <- length(sel)
            if (nGenesDn[iGeneSet] > 0) {
              gsStatsDn[iGeneSet] <- calcGeneSetStat(sel,
                                                     method)
            }
          }
        }
      }
      else if (method %in% c("stouffer", "reporter", "tailStrength")) {
        indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
        statsGenesInSet <- statsContrast[indGenesInSet]
        nGenes[iGeneSet] <- length(statsGenesInSet)
        if (nGenes[iGeneSet] > 0) {
          if (exists("signsContrast")) {
            nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet] > 
                                        0)
            nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet] < 
                                        0)
          }
          if (statType == "p-signed") {
            statsGenesInSetTestUp <- statsContrastTestUp[indGenesInSet]
            statsGenesInSetTestDn <- statsContrastTestDn[indGenesInSet]
            statsGenesInSetTestUp <- statsGenesInSetTestUp[!is.na(statsGenesInSetTestUp)]
            statsGenesInSetTestDn <- statsGenesInSetTestDn[!is.na(statsGenesInSetTestDn)]
            gsStatsAllTestUp[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestUp,
                                                          method)
            gsStatsAllTestDn[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestDn,
                                                          method)
          }
          gsStatsAbs[iGeneSet] <- calcGeneSetStat(abs(statsGenesInSet),
                                                  method)
          if (statType == "p-signed" & signMethod != 
              "sampleperm") {
            signsGenesInSet <- signsContrast[indGenesInSet]
            sel <- statsGenesInSet[signsGenesInSet > 
                                     0]
            nGenesUp[iGeneSet] <- length(sel)
            if (nGenesUp[iGeneSet] > 0) {
              gsStatsUp[iGeneSet] <- calcGeneSetStat(sel,
                                                     method)
            }
            sel <- abs(statsGenesInSet[signsGenesInSet < 
                                         0])
            nGenesDn[iGeneSet] <- length(sel)
            if (nGenesDn[iGeneSet] > 0) {
              gsStatsDn[iGeneSet] <- calcGeneSetStat(sel,
                                                     method)
            }
          }
        }
      }
      else if (method == "wilcoxon") {
        indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
        statsGenesInSet <- statsContrast[indGenesInSet]
        statsNotInSet <- statsContrast[!indGenesInSet]
        if (statType %in% c("p-signed", "F-signed")) {
          signsGenesInSet <- signsContrast[indGenesInSet]
          signsNotInSet <- signsContrast[!indGenesInSet]
        }
        nGenes[iGeneSet] <- length(statsGenesInSet)
        if (nGenes[iGeneSet] > 0) {
          if (exists("signsContrast")) {
            nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet] > 
                                        0)
            nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet] < 
                                        0)
          }
          if (statType == "p-signed") {
            statsGenesInSetTestUp <- statsContrastTestUp[indGenesInSet]
            statsGenesInSetTestDn <- statsContrastTestDn[indGenesInSet]
            statsGenesInSetTestUp <- statsGenesInSetTestUp[!is.na(statsGenesInSetTestUp)]
            statsGenesInSetTestDn <- statsGenesInSetTestDn[!is.na(statsGenesInSetTestDn)]
            statsNotInSetTestUp <- statsContrastTestUp[!indGenesInSet]
            statsNotInSetTestDn <- statsContrastTestDn[!indGenesInSet]
            statsNotInSetTestUp <- statsNotInSetTestUp[!is.na(statsNotInSetTestUp)]
            statsNotInSetTestDn <- statsNotInSetTestDn[!is.na(statsNotInSetTestDn)]
            tmp <- calcGeneSetStat(statsGenesInSetTestUp,
                                   "wilcoxon_less", statsNotInSetTestUp)
            gsStatsAllTestUp[iGeneSet] <- tmp[1]
            pValuesAllUp[iGeneSet] <- tmp[2]
            tmp <- calcGeneSetStat(statsGenesInSetTestDn,
                                   "wilcoxon_less", statsNotInSetTestDn)
            gsStatsAllTestDn[iGeneSet] <- tmp[1]
            pValuesAllDn[iGeneSet] <- tmp[2]
          }
          else if (statType == "t") {
            tmp <- calcGeneSetStat(statsGenesInSet, "wilcoxon_two.sided",
                                   statsNotInSet)
            gsStatsAll[iGeneSet] <- tmp[1]
            tmp <- calcGeneSetStat(statsGenesInSet, "wilcoxon_greater",
                                   statsNotInSet)
            pValuesAllUp[iGeneSet] <- tmp[2]
            tmp <- calcGeneSetStat(statsGenesInSet, "wilcoxon_less",
                                   statsNotInSet)
            pValuesAllDn[iGeneSet] <- tmp[2]
          }
          if (statType %in% c("p", "p-signed")) {
            tmp <- calcGeneSetStat(abs(statsGenesInSet),
                                   "wilcoxon_less", abs(statsNotInSet))
            gsStatsAbs[iGeneSet] <- tmp[1]
            pValuesAbs[iGeneSet] <- tmp[2]
          }
          else if (statType %in% c("t", "F", "F-signed")) {
            tmp <- calcGeneSetStat(abs(statsGenesInSet),
                                   "wilcoxon_greater", abs(statsNotInSet))
            gsStatsAbs[iGeneSet] <- tmp[1]
            pValuesAbs[iGeneSet] <- tmp[2]
          }
          if (statType %in% c("t", "p-signed", "F-signed") & 
              signMethod != "sampleperm") {
            if (statType == "t") {
              sel <- statsGenesInSet[statsGenesInSet > 
                                       0]
            }
            else {
              sel <- statsGenesInSet[signsGenesInSet > 
                                       0]
            }
            nGenesUp[iGeneSet] <- length(sel)
            if (nGenesUp[iGeneSet] > 0) {
              if (statType == "p-signed") {
                tmp <- calcGeneSetStat(sel, "wilcoxon_less",
                                       statsNotInSet[signsNotInSet > 0])
                gsStatsUp[iGeneSet] <- tmp[1]
                pValuesUp[iGeneSet] <- tmp[2]
              }
              else if (statType == "t") {
                tmp <- calcGeneSetStat(sel, "wilcoxon_greater",
                                       statsNotInSet[statsNotInSet > 0])
                gsStatsUp[iGeneSet] <- tmp[1]
                pValuesUp[iGeneSet] <- tmp[2]
              }
              else if (statType == "F-signed") {
                tmp <- calcGeneSetStat(sel, "wilcoxon_greater",
                                       statsNotInSet[signsNotInSet > 0])
                gsStatsUp[iGeneSet] <- tmp[1]
                pValuesUp[iGeneSet] <- tmp[2]
              }
            }
            if (statType == "t") {
              sel <- abs(statsGenesInSet[statsGenesInSet < 
                                           0])
            }
            else {
              sel <- abs(statsGenesInSet[signsGenesInSet < 
                                           0])
            }
            nGenesDn[iGeneSet] <- length(sel)
            if (nGenesDn[iGeneSet] > 0) {
              if (statType == "p-signed") {
                tmp <- calcGeneSetStat(sel, "wilcoxon_less", 
                                       abs(statsNotInSet[signsNotInSet < 0]))
                gsStatsDn[iGeneSet] <- tmp[1]
                pValuesDn[iGeneSet] <- tmp[2]
              }
              else if (statType == "t") {
                tmp <- calcGeneSetStat(sel, "wilcoxon_greater", 
                                       abs(statsNotInSet[statsNotInSet < 0]))
                gsStatsDn[iGeneSet] <- tmp[1]
                pValuesDn[iGeneSet] <- tmp[2]
              }
              else if (statType == "F-signed") {
                tmp <- calcGeneSetStat(sel, "wilcoxon_greater", 
                                       abs(statsNotInSet[signsNotInSet < 0]))
                gsStatsDn[iGeneSet] <- tmp[1]
                pValuesDn[iGeneSet] <- tmp[2]
              }
            }
          }
        }
      }
      else if (method == "wilcoxon_fast") {
        indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
        statsGenesInSet <- statsContrast[indGenesInSet]
        statsNotInSet <- statsContrast[!indGenesInSet]
        if (statType %in% c("p-signed", "F-signed")) {
          signsGenesInSet <- signsContrast[indGenesInSet]
          signsNotInSet <- signsContrast[!indGenesInSet]
        }
        nGenes[iGeneSet] <- length(statsGenesInSet)
        if (nGenes[iGeneSet] > 0) {
          if (exists("signsContrast")) {
            nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet] > 
                                        0)
            nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet] < 
                                        0)
          }
          if (statType == "p-signed") {
            statsGenesInSetTestUp <- statsContrastTestUp[indGenesInSet]
            statsGenesInSetTestDn <- statsContrastTestDn[indGenesInSet]
            statsGenesInSetTestUp <- statsGenesInSetTestUp[!is.na(statsGenesInSetTestUp)]
            statsGenesInSetTestDn <- statsGenesInSetTestDn[!is.na(statsGenesInSetTestDn)]
            statsNotInSetTestUp <- statsContrastTestUp[!indGenesInSet]
            statsNotInSetTestDn <- statsContrastTestDn[!indGenesInSet]
            statsNotInSetTestUp <- statsNotInSetTestUp[!is.na(statsNotInSetTestUp)]
            statsNotInSetTestDn <- statsNotInSetTestDn[!is.na(statsNotInSetTestDn)]
            gsStatsAllTestUp[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestUp,
                                                          "wilcoxon_fast", statsNotInSetTestUp)
            gsStatsAllTestDn[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestDn,
                                                          "wilcoxon_fast", statsNotInSetTestDn)
          }
          else if (statType == "t") {
            gsStatsAll[iGeneSet] <- calcGeneSetStat(statsGenesInSet,
                                                    "wilcoxon_fast", statsNotInSet)
          }
          if (statType %in% c("p", "p-signed")) {
            gsStatsAbs[iGeneSet] <- calcGeneSetStat(abs(statsGenesInSet),
                                                    "wilcoxon_fast", abs(statsNotInSet))
          }
          else if (statType %in% c("t", "F", "F-signed")) {
            gsStatsAbs[iGeneSet] <- calcGeneSetStat(abs(statsGenesInSet),
                                                    "wilcoxon_fast", abs(statsNotInSet))
          }
          if (statType %in% c("t", "p-signed", "F-signed") & 
              signMethod != "sampleperm") {
            if (statType == "t") {
              sel <- statsGenesInSet[statsGenesInSet > 
                                       0]
            }
            else {
              sel <- statsGenesInSet[signsGenesInSet > 
                                       0]
            }
            nGenesUp[iGeneSet] <- length(sel)
            if (nGenesUp[iGeneSet] > 0) {
              if (statType == "p-signed") {
                gsStatsUp[iGeneSet] <- calcGeneSetStat(sel,
                                                       "wilcoxon_fast", statsNotInSet[signsNotInSet > 
                                                                                        0])
              }
              else if (statType == "t") {
                gsStatsUp[iGeneSet] <- calcGeneSetStat(sel,
                                                       "wilcoxon_fast", statsNotInSet[statsNotInSet > 
                                                                                        0])
              }
              else if (statType == "F-signed") {
                gsStatsUp[iGeneSet] <- calcGeneSetStat(sel,
                                                       "wilcoxon_fast", statsNotInSet[signsNotInSet > 
                                                                                        0])
              }
            }
            if (statType == "t") {
              sel <- abs(statsGenesInSet[statsGenesInSet < 
                                           0])
            }
            else {
              sel <- abs(statsGenesInSet[signsGenesInSet < 
                                           0])
            }
            nGenesDn[iGeneSet] <- length(sel)
            if (nGenesDn[iGeneSet] > 0) {
              if (statType == "p-signed") {
                gsStatsDn[iGeneSet] <- calcGeneSetStat(sel,
                                                       "wilcoxon_fast", abs(statsNotInSet[signsNotInSet < 
                                                                                            0]))
              }
              else if (statType == "t") {
                gsStatsDn[iGeneSet] <- calcGeneSetStat(sel,
                                                       "wilcoxon_fast", abs(statsNotInSet[statsNotInSet < 
                                                                                            0]))
              }
              else if (statType == "F-signed") {
                gsStatsDn[iGeneSet] <- calcGeneSetStat(sel,
                                                       "wilcoxon_fast", abs(statsNotInSet[signsNotInSet < 
                                                                                            0]))
              }
            }
          }
        }
      }
      else if (method %in% c("mean", "median", "sum")) {
        indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
        statsGenesInSet <- statsContrast[indGenesInSet]
        if (statType %in% c("p-signed", "F-signed")) {
          signsGenesInSet <- signsContrast[indGenesInSet]
        }
        nGenes[iGeneSet] <- length(statsGenesInSet)
        if (nGenes[iGeneSet] > 0) {
          if (exists("signsContrast")) {
            nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet] > 
                                        0)
            nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet] < 
                                        0)
          }
          if (statType == "p-signed") {
            statsGenesInSetTestUp <- statsContrastTestUp[indGenesInSet]
            statsGenesInSetTestDn <- statsContrastTestDn[indGenesInSet]
            statsGenesInSetTestUp <- statsGenesInSetTestUp[!is.na(statsGenesInSetTestUp)]
            statsGenesInSetTestDn <- statsGenesInSetTestDn[!is.na(statsGenesInSetTestDn)]
            gsStatsAllTestUp[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestUp,
                                                          method)
            gsStatsAllTestDn[iGeneSet] <- calcGeneSetStat(statsGenesInSetTestDn,
                                                          method)
          }
          else if (statType == "t") {
            gsStatsAll[iGeneSet] <- calcGeneSetStat(statsGenesInSet, 
                                                    method)
          }
          gsStatsAbs[iGeneSet] <- calcGeneSetStat(abs(statsGenesInSet), 
                                                  method)
          if (statType %in% c("p-signed", "t", "F-signed") & 
              signMethod != "sampleperm") {
            if (statType == "t") {
              sel <- statsGenesInSet[statsGenesInSet > 
                                       0]
            }
            else {
              sel <- statsGenesInSet[signsGenesInSet > 
                                       0]
            }
            nGenesUp[iGeneSet] <- length(sel)
            if (nGenesUp[iGeneSet] > 0) {
              gsStatsUp[iGeneSet] <- calcGeneSetStat(sel, 
                                                     method)
            }
            if (statType == "t") {
              sel <- abs(statsGenesInSet[statsGenesInSet < 
                                           0])
            }
            else {
              sel <- abs(statsGenesInSet[signsGenesInSet < 
                                           0])
            }
            nGenesDn[iGeneSet] <- length(sel)
            if (nGenesDn[iGeneSet] > 0) {
              gsStatsDn[iGeneSet] <- calcGeneSetStat(sel, 
                                                     method)
            }
          }
        }
      }
      else if (method == "maxmean") {
        indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
        statsGenesInSet <- statsContrast[indGenesInSet]
        nGenes[iGeneSet] <- length(statsGenesInSet)
        if (nGenes[iGeneSet] > 0) {
          if (exists("signsContrast")) {
            nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet] > 
                                        0)
            nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet] < 
                                        0)
          }
          gsStatsAbs[iGeneSet] <- calcGeneSetStat(statsGenesInSet, 
                                                  method)
        }
      }
      else if (method == "gsea") {
        statsContrastSorted <- sort(statsContrast, decreasing = TRUE)
        indGenesInSet <- which(names(statsContrastSorted) %in% 
                                 gsc[[iGeneSet]])
        nGenes[iGeneSet] <- length(indGenesInSet)
        if (nGenes[iGeneSet] > 0) {
          if (exists("signsContrast")) {
            nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet] > 
                                        0)
            nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet] < 
                                        0)
          }
          gsStatsAll[iGeneSet] <- calcGeneSetStat(indGenesInSet, 
                                                  method, statsContrastSorted, gseaParam)
        }
      }
      else if (method == "page") {
        indGenesInSet <- names(statsContrast) %in% gsc[[iGeneSet]]
        statsGenesInSet <- statsContrast[indGenesInSet]
        nGenes[iGeneSet] <- length(statsGenesInSet)
        if (nGenes[iGeneSet] > 0) {
          if (exists("signsContrast")) {
            nGenesUp[iGeneSet] <- sum(signsContrast[indGenesInSet] > 
                                        0)
            nGenesDn[iGeneSet] <- sum(signsContrast[indGenesInSet] < 
                                        0)
          }
          gsStatsAll[iGeneSet] <- calcGeneSetStat(statsGenesInSet, 
                                                  method, statsContrast)
        }
      }
    }
    res$statsAll <- cbind(res$statsAll, gsStatsAll)
    res$statsAllTestUp <- cbind(res$statsAllTestUp, gsStatsAllTestUp)
    res$statsAllTestDn <- cbind(res$statsAllTestDn, gsStatsAllTestDn)
    res$statsAbs <- cbind(res$statsAbs, gsStatsAbs)
    res$statsUp <- cbind(res$statsUp, gsStatsUp)
    res$statsDn <- cbind(res$statsDn, gsStatsDn)
    res$nGenes <- cbind(res$nGenes, nGenes)
    res$nGenesUp <- cbind(res$nGenesUp, nGenesUp)
    res$nGenesDn <- cbind(res$nGenesDn, nGenesDn)
    res$pValuesAll <- cbind(res$pValuesAll, pValuesAll)
    res$pValuesAllUp <- cbind(res$pValuesAllUp, pValuesAllUp)
    res$pValuesAllDn <- cbind(res$pValuesAllDn, pValuesAllDn)
    res$pValuesAbs <- cbind(res$pValuesAbs, pValuesAbs)
    res$pValuesUp <- cbind(res$pValuesUp, pValuesUp)
    res$pValuesDn <- cbind(res$pValuesDn, pValuesDn)
  }
  if (method == "fisher") 
    res$statName = "X2"
  else if (method == "stouffer") 
    res$statName = "Z"
  else if (method == "tailStrength") 
    res$statName = "TS"
  else if (method == "reporter") 
    res$statName = "Z"
  else if (method == "wilcoxon") 
    res$statName = "W"
  else if (method == "wilcoxon_fast") 
    res$statName = "W"
  else if (method == "mean") 
    res$statName = "mean"
  else if (method == "median") 
    res$statName = "median"
  else if (method == "sum") 
    res$statName = "sum"
  else if (method == "maxmean") 
    res$statName = "maxmean"
  else if (method == "gsea") 
    res$statName = "ES"
  else if (method == "page") 
    res$statName = "Z"
  return(res)
}