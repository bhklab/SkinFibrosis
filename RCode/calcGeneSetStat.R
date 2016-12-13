

## function adapted from Piano package: https://bioconductor.org/packages/release/bioc/html/piano.html

calcGeneSetStat <- function (selectedStats, method, statistics = NULL, gseaParam)
{
        if (method == "fisher") {
            geneSetStatistic <- 2 * (sum(-1 * log(selectedStats)))
        }
        else if (method %in% c("stouffer", "reporter")) {
            geneSetStatistic <- sum(qnorm(selectedStats, lower.tail = FALSE))/sqrt(length(selectedStats))
        }
        else if (method == "tailStrength") {
            m <- length(selectedStats)
            geneSetStatistic <- (1/m) * sum(1 - sort(selectedStats) *
            (m + 1)/(1:m))
        }
        else if (method == "wilcoxon_less") {
            geneSetStatistic <- unlist(wilcox.test(selectedStats,
            statistics, alternative = "less", conf.int = FALSE)[c(1,
            3)])
        }
        else if (method == "wilcoxon_greater") {
            geneSetStatistic <- unlist(wilcox.test(selectedStats,
            statistics, alternative = "greater", conf.int = FALSE)[c(1,
            3)])
        }
        else if (method == "wilcoxon_two.sided") {
            geneSetStatistic <- unlist(wilcox.test(selectedStats,
            statistics, alternative = "two.sided", conf.int = FALSE)[c(1,
            3)])
        }
        else if (method == "wilcoxon_fast") {
            m <- length(selectedStats)
            geneSetStatistic <- sum(rank(c(selectedStats, statistics))[1:m]) -
            m * (m + 1)/2
        }
        else if (method == "mean") {
            geneSetStatistic <- mean(selectedStats)
        }
        else if (method == "median") {
            geneSetStatistic <- median(selectedStats)
        }
        else if (method == "sum") {
            geneSetStatistic <- sum(selectedStats)
        }
        else if (method == "maxmean") {
            m <- length(selectedStats)
            sPlus <- sum(selectedStats[selectedStats > 0])/m
            sMinus <- -sum(selectedStats[selectedStats < 0])/m
            geneSetStatistic <- max(c(sPlus, sMinus))
        }
        else if (method == "gsea") {
            S <- selectedStats
            r <- statistics
            p <- gseaParam
            m <- length(S)
            N <- length(r)
            NR <- (sum(abs(r[S])^p))
            Pprev <- 0
            P <- rep(NA, N)
            for (i in 1:N) {
                if (i %in% S) {
                    if (r[i] == 0) {
                        P[i] <- Pprev + 0
                    }
                    else {
                        P[i] <- Pprev + abs(r[i])^p/NR
                    }
                }
                else {
                    P[i] <- Pprev - 1/(N - m)
                }
                Pprev <- P[i]
            }
            if (max(P) > -min(P)) {
                geneSetStatistic <- max(P)
            }
            else {
                geneSetStatistic <- min(P)
            }
        }
        else if (method == "page") {
            mu <- mean(statistics)
            delta <- sd(statistics)
            Sm <- mean(selectedStats)
            m <- length(selectedStats)
            geneSetStatistic <- (Sm - mu) * sqrt(m)/delta
        }
        return(geneSetStatistic)
    }