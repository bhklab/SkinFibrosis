## function adapted from Piano package: https://bioconductor.org/packages/release/bioc/html/piano.html


GSCsignificanceBatch <- function (statistics, statType, signs, gsc, statMethod, signMethod,
          permStatistics, permSigns, nGenes, nGenesUp, nGenesDn, gsStatsAll, 
          gsStatsAllTestUp, gsStatsAllTestDn, gsStatsAbs, gsStatsUp, 
          gsStatsDn, nPerm, gseaParam, ncpus) 
{
  if (signMethod == "geneperm" | (statMethod == "reporter" & 
                                  signMethod == "distribution")) {
    if (ncpus > 1) {
      
      
      if (!requireNamespace("snowfall", quietly = TRUE)) 
        stop("package snowfall is missing")
      tmp <- capture.output(suppressMessages(snowfall::sfInit(parallel = T, 
                                                              cpus = ncpus)))
      tmp <- snowfall::sfLapply(1:ncpus, GSCstatGenePerm, 
                                statistics, signs, gsc, statType, statMethod, 
                                nGenes, nGenesUp, nGenesDn, nPerm/ncpus, gseaParam)
      suppressMessages(snowfall::sfStop())
      gsStatsAllPerm <- vector()
      gsStatsAllTestUpPerm <- vector()
      gsStatsAllTestDnPerm <- vector()
      gsStatsAbsPerm <- vector()
      gsStatsUpPerm <- vector()
      gsStatsDnPerm <- vector()
      for (i in 1:ncpus) {
        gsStatsAllPerm <- cbind(gsStatsAllPerm, tmp[[i]]$gsStatsAllPerm[[1]])
        gsStatsAllTestUpPerm <- cbind(gsStatsAllTestUpPerm, 
                                      tmp[[i]]$gsStatsAllTestUpPerm[[1]])
        gsStatsAllTestDnPerm <- cbind(gsStatsAllTestDnPerm, 
                                      tmp[[i]]$gsStatsAllTestDnPerm[[1]])
        gsStatsAbsPerm <- cbind(gsStatsAbsPerm, tmp[[i]]$gsStatsAbsPerm[[1]])
        gsStatsUpPerm <- cbind(gsStatsUpPerm, tmp[[i]]$gsStatsUpPerm[[1]])
        gsStatsDnPerm <- cbind(gsStatsDnPerm, tmp[[i]]$gsStatsDnPerm[[1]])
      }
      gsStatsAllPerm <- list(gsStatsAllPerm)
      gsStatsAllTestUpPerm <- list(gsStatsAllTestUpPerm)
      gsStatsAllTestDnPerm <- list(gsStatsAllTestDnPerm)
      gsStatsAbsPerm <- list(gsStatsAbsPerm)
      gsStatsUpPerm <- list(gsStatsUpPerm)
      gsStatsDnPerm <- list(gsStatsDnPerm)
    }
    else {
      res <- GSCstatGenePerm(NULL, statistics, signs, gsc, 
                             statType, statMethod, nGenes, nGenesUp, nGenesDn, 
                             nPerm, gseaParam)
      gsStatsAllPerm <- res$gsStatsAllPerm
      gsStatsAllTestUpPerm <- res$gsStatsAllTestUpPerm
      gsStatsAllTestDnPerm <- res$gsStatsAllTestDnPerm
      gsStatsAbsPerm <- res$gsStatsAbsPerm
      gsStatsUpPerm <- res$gsStatsUpPerm
      gsStatsDnPerm <- res$gsStatsDnPerm
    }
  }
  else if (signMethod == "sampleperm") {
    res <- GSCstatSamplePerm(permStatistics, permSigns, statType, 
                             gsc, statMethod, gseaParam)
    gsStatsAllPerm <- res$gsStatsAllPerm
    gsStatsAllTestUpPerm <- res$gsStatsAllTestUpPerm
    gsStatsAllTestDnPerm <- res$gsStatsAllTestDnPerm
    gsStatsAbsPerm <- res$gsStatsAbsPerm
    gsStatsUpPerm <- res$gsStatsUpPerm
    gsStatsDnPerm <- res$gsStatsDnPerm
  }
  else {
    gsStatsAllPerm <- NA
    gsStatsAllTestUpPerm <- NA
    gsStatsAllTestDnPerm <- NA
    gsStatsAbsPerm <- NA
    gsStatsUpPerm <- NA
    gsStatsDnPerm <- NA
  }
  if (signMethod == "geneperm") {
    res <- pvalFromFractionGenePerm(statType, statMethod, 
                                    nGenes, nGenesUp, nGenesDn, gsStatsAll, gsStatsAllTestUp, 
                                    gsStatsAllTestDn, gsStatsAbs, gsStatsUp, gsStatsDn, 
                                    gsStatsAllPerm, gsStatsAllTestUpPerm, gsStatsAllTestDnPerm, 
                                    gsStatsAbsPerm, gsStatsUpPerm, gsStatsDnPerm)
  }
  else if (signMethod == "sampleperm") {
    res <- pvalFromFractionSamplePerm(statType, statMethod, 
                                      nGenes, nGenesUp, nGenesDn, gsStatsAll, gsStatsAllTestUp, 
                                      gsStatsAllTestDn, gsStatsAbs, gsStatsUp, gsStatsDn, 
                                      gsStatsAllPerm, gsStatsAllTestUpPerm, gsStatsAllTestDnPerm, 
                                      gsStatsAbsPerm, gsStatsUpPerm, gsStatsDnPerm)
  }
  if (signMethod == "distribution") {
    res <- pvalFromDistribution(statType, statMethod, nGenes, 
                                nGenesUp, nGenesDn, gsStatsAll, gsStatsAllTestUp, 
                                gsStatsAllTestDn, gsStatsAbs, gsStatsUp, gsStatsDn, 
                                gsStatsAllPerm, gsStatsAllTestUpPerm, gsStatsAllTestDnPerm, 
                                gsStatsAbsPerm, gsStatsUpPerm, gsStatsDnPerm)
  }
  res$gsStatsAllPerm <- gsStatsAllPerm
  return(res)
}