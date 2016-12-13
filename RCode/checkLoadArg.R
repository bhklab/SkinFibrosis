
## function adapted from Piano package: https://bioconductor.org/packages/release/bioc/html/piano.html


checkLoadArg <- function (statistics, signs, statMethod, signMethod, adjMethod,
          gsc, gsSize, permStatistics, permSigns, nPerm, gseaParam, 
          ncpus, verbose) 
{
  if (class(verbose) != "logical") 
    stop("argument verbose has to be TRUE or FALSE")
  if (class(ncpus) != "numeric") 
    stop("argument ncpus should be an integer")
  if (length(ncpus) != 1) 
    stop("argument ncpus should be an integer")
  tmp <- try(statMethod <- match.arg(statMethod, c("fisher", 
                                                   "stouffer", "reporter", "tailStrength", "wilcoxon", "mean", 
                                                   "median", "sum", "maxmean", "gsea", "page"), several.ok = FALSE), 
             silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument geneSetStat set to unknown method")
  }
  tmp <- try(signMethod <- match.arg(signMethod, c("geneSampling", 
                                                   "samplePermutation", "nullDist"), several.ok = FALSE), 
             silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument signifMethod set to unknown method")
  }
  if (signMethod == "geneSampling") 
    signMethod <- "geneperm"
  if (signMethod == "samplePermutation") 
    signMethod <- "sampleperm"
  if (signMethod == "nullDist") 
    signMethod <- "distribution"
  if (signMethod == "sampleperm" & is.null(permStatistics)) {
    stop("signifMethod='samplePermutation' but no permStats given, can not perform sample permutation")
  }
  if (signMethod == "geneperm" & !is.null(permStatistics)) {
    warning("signifMethod='geneSampling' and argument permStats given, will not use permStats")
  }
  if (signMethod == "distribution" & !statMethod %in% c("fisher", 
                                                        "stouffer", "reporter", "wilcoxon", "page")) {
    stop(paste("signifMethod='nullDist' is not allowed for geneSetStat='", 
               statMethod, "'", sep = ""))
  }
  tmp <- try(adjMethod <- match.arg(adjMethod, c("holm", "hochberg", 
                                                 "hommel", "bonferroni", "BH", "BY", "fdr", "none"), several.ok = FALSE), 
             silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument adjMethod set to unknown method")
  }
  if (!adjMethod %in% c("fdr", "none") & statMethod == "gsea") {
    statMethod <- "fdr"
    warning("adjMethod can only be 'fdr' for statMethod='gsea', using fdr despite different user setting")
  }
  if (length(nPerm) != 1) 
    stop("length of argument nPerm has to be 1")
  if (nPerm < 100) 
    stop("argument nPerm has to be >100")
  if (nPerm%%ncpus != 0) 
    stop("argument ncpus should be set so there is an integer x such that x*ncpus=nPerm")
  if (gseaParam < 1) 
    stop("gseaParam has to be larger than 0")
  tmp <- try(statistics <- as.matrix(statistics), silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument geneLevelStats could not be converted into a matrix")
  }
  if (ncol(statistics) != 1) 
    stop("geneLevelStats should be a vector or only contain one column")
  if (sum(is.na(statistics)) > 0) {
    stop("NA values not allowed in geneLevelStats")
  }
  if (min(statistics) >= 0 & max(statistics) <= 1) {
    statType <- "p"
  }
  else if (sign(min(statistics)) != sign(max(statistics))) {
    statType <- "t"
  }
  else {
    statType <- "F"
  }
  if (statType == "p") {
    statistics[statistics < 1e-100] <- 1e-100
    statistics[statistics > 1] <- 1
  }
  if (statMethod %in% c("fisher", "stouffer", "reporter", "tailStrength") & 
      statType != "p") {
    stop(paste("geneLevelStats does not lie in [0,1], geneLevelStats can only be p-values for geneSetStat='", 
               statMethod, "'", sep = ""))
  }
  if (statMethod %in% c("maxmean", "gsea", "page") & statType != 
      "t") {
    stop(paste("geneLevelStats has to contain both positive and negative scores for geneSetStat='", 
               statMethod, "'", sep = ""))
  }
  if (statMethod == "page") {
    for (i in 1:ncol(statistics)) {
      if (sd(statistics[, i]) > 1e+300) {
        stop("standard deviation of geneLevelStats is close to infinity, for geneSetStat 'page' all gene-set statistics will be zero, change geneLevelStats or geneSetStat")
      }
    }
  }
  tmp <- rownames(statistics)
  if (length(tmp) > length(unique(tmp))) 
    warning("found duplicates in rownames(geneLevelStats), all values will be used for calculation of gene set statistics")
  info <- list()
  info$nGenesStatistics <- length(rownames(statistics))
  info$nContrasts <- ncol(statistics)
  if (!is.null(signs) & statType != "t") {
    tmp <- try(signs <- as.matrix(signs), silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument directions could not be converted into a matrix")
    }
    if (ncol(signs) != 1) 
      stop("argument directions should be a vector or only contain one column")
    if (signMethod == "sampleperm" & is.null(permSigns)) {
      stop("for signifMethod='samplePermutation', both directions and permDirections have to be given, or vice versa")
    }
    if (sum(is.na(signs)) > 0) {
      stop("NA values not allowed in directions")
    }
    if (nrow(statistics) != nrow(signs) | ncol(statistics) != 
        ncol(signs)) {
      stop("dimensions of geneLevelStats and directions are not the same")
    }
    if (is.null(rownames(statistics)) | is.null(rownames(signs))) {
      stop("rownames of geneLevelStats and directions do not match")
    }
    if (!identical(rownames(statistics), rownames(signs))) {
      stop("rownames of geneLevelStats and directions do not match")
    }
    signs <- sign(signs)
    if (statType == "p") 
      statType <- "p-signed"
    if (statType == "F") 
      statType <- "F-signed"
  }
  else if (!is.null(signs) & statType == "t") {
    warning("geneLevelStats are t-like and do not require information given by argument directions, argument directions will not be used")
    signs <- "none"
  }
  else {
    signs <- "none"
  }
  if (class(gsc) != "GSC") 
    stop("argument gsc should have class 'GSC' as output from loadGSC()")
  addInfo <- gsc$addInfo
  gsc <- gsc$gsc
  tmp <- vector()
  for (iGeneSet in 1:length(gsc)) {
    tmp <- c(tmp, gsc[[iGeneSet]][!gsc[[iGeneSet]] %in% rownames(statistics)])
    gsc[[iGeneSet]] <- gsc[[iGeneSet]][gsc[[iGeneSet]] %in% 
                                         rownames(statistics)]
  }
  info$removedGenesGSC <- length(unique(tmp))
  tmp <- length(gsc)
  gsc <- gsc[unlist(lapply(gsc, length)) > 0]
  info$removedGSnoGenes <- tmp - length(gsc)
  if (length(gsSize) != 2) 
    stop("length(gsSizeLim) has to equal 2")
  if (gsSize[1] > gsSize[2]) 
    stop("gsSizeLim[1] has to be smaller or equal to gsSizeLim[2]")
  if (gsSize[1] < 1) 
    stop("gsSizeLim[1] has to be larger than 0")
  tmp <- length(gsc)
  gsc <- gsc[unlist(lapply(gsc, length)) >= gsSize[1] & unlist(lapply(gsc, 
                                                                      length)) <= gsSize[2]]
  info$removedGSsizeLimit <- tmp - length(gsc)
  if (class(addInfo) != "character") {
    tmp <- sum(addInfo[, 1] %in% names(gsc))
    info$nGeneSetsWithAddInfo <- tmp
  }
  else {
    info$nGeneSetsWithAddInfo <- 0
  }
  info$nGenesGSC <- length(unique(unlist(gsc)))
  info$nGeneSets <- length(gsc)
  if (!is.null(permStatistics) & signMethod == "sampleperm") {
    permStatistics <- list(permStatistics)
    permSigns <- list(permSigns)
    if (!is.null(permSigns) & statType == "t") {
      warning("geneLevelStats are t-like and do not require information given by argument permDirections, argument permDirections will not be used")
    }
    if (is.null(signs) & !is.null(permSigns)) {
      if (statType != "t") 
        stop("for signifMethod='samplePermutation', both directions and permDirections have to be given, or vice versa")
    }
    for (i in 1:length(permStatistics)) {
      tmp <- try(permStatistics[[i]] <- as.matrix(permStatistics[[i]]), 
                 silent = TRUE)
      if (class(tmp) == "try-error") {
        stop("permStats could not be converted into a matrix")
      }
      if (sum(is.na(permStatistics[[i]])) > 0) {
        stop("NA values not allowed in permStats")
      }
      if (min(permStatistics[[i]]) >= 0 & max(permStatistics[[i]]) <= 
          1) {
        tmp <- "p"
      }
      else if (sign(min(permStatistics[[i]])) != sign(max(permStatistics[[i]]))) {
        tmp <- "t"
      }
      else {
        tmp <- "F"
      }
      if (tmp == "p" & !statType %in% c("p", "p-signed")) 
        stop("geneLevelStats and permStats has to contain the same type of statistics")
      if (tmp == "F" & !statType %in% c("F", "F-signed")) 
        stop("geneLevelStats and permStats has to contain the same type of statistics")
      if (tmp == "t" & !statType == "t") 
        stop("geneLevelStats and permStats has to contain the same type of statistics")
      if (tmp == "p") {
        permStatistics[[i]][permStatistics[[i]] < 1e-100] <- 1e-100
        permStatistics[[i]][permStatistics[[i]] > 1] <- 1
      }
      if (statMethod %in% c("fisher", "stouffer", "reporter", 
                            "tailStrength") & tmp != "p") {
        stop(paste("permStats does not lie in [0,1], permStats can only be p-values for geneSetStat='", 
                   statMethod, "'", sep = ""))
      }
      if (statMethod %in% c("maxmean", "gsea", "page") & 
          tmp != "t") {
        stop(paste("permStats has to contain both positive and negative scores for geneSetStat='", 
                   statMethod, "'", sep = ""))
      }
      if (nrow(statistics) != nrow(permStatistics[[i]])) {
        stop("permStats does not have the same number of rows as geneLevelStats")
      }
      if (!identical(rownames(statistics), rownames(permStatistics[[i]]))) {
        stop("rownames of geneLevelStats and permStats do not match")
      }
      if (!is.null(permSigns) & statType != "t") {
        tmp <- try(permSigns[[i]] <- as.matrix(permSigns[[i]]), 
                   silent = TRUE)
        if (class(tmp) == "try-error") {
          stop("permDirections could not be converted into a matrix")
        }
        if (sum(is.na(signs[[i]])) > 0) {
          stop("NA values not allowed in permDirections")
        }
        if (nrow(permStatistics[[i]]) != nrow(permSigns[[i]]) | 
            ncol(permStatistics[[i]]) != ncol(permSigns[[i]])) {
          stop("dimensions of permStats and permDirections are not the same")
        }
        if (!identical(rownames(permStatistics[[i]]), 
                       rownames(permSigns[[i]]))) {
          stop("rownames of permStats and permDirections do not match")
        }
        permSigns[[i]] <- sign(permSigns[[i]])
      }
      else {
        permSigns <- "none"
      }
    }
    info$nSamplePermutations <- ncol(permStatistics[[1]])
    nPerm <- ncol(permStatistics[[1]])
  }
  else {
    permStatistics <- "none"
  }
  res <- list()
  res$statistics <- statistics
  res$statType <- statType
  res$statMethod <- statMethod
  res$signMethod <- signMethod
  res$adjMethod <- adjMethod
  res$gsc <- gsc
  res$permStatistics <- permStatistics
  res$addInfo <- addInfo
  res$info <- info
  res$dataSubsets <- "Not used"
  res$nPerm <- nPerm
  res$signs <- signs
  res$permSigns <- permSigns
  return(res)
}