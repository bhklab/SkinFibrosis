

querySigKEGG <- function() {
    
sourceName="KEGG"
aa <- read.csv("../Data/CPDB_pathways_genes-human.tab", sep="\t", comment.char="#")
aa <- aa[aa$source==sourceName,] 
aa$pathway <- as.character(aa$pathway)
aa$hgnc_symbol_ids <- as.character(aa$hgnc_symbol_ids)
aa$hgnc_symbol_ids <- toupper(aa$hgnc_symbol_ids)


## ppar
ssup <- unlist(strsplit(as.character(aa[iup,]$hgnc_symbol_ids), ","))
##glyco
ssdn <- unlist(strsplit(as.character(aa[idn,]$hgnc_symbol_ids), ","))
length(intersect(ssup, ssdn))
#[1] 2

xx <- intersect(ssup, ssdn)
#[1] "PCK2" "PCK1"
ssup <- setdiff(ssup, xx) ## down 
ssdn <- setdiff(ssdn,xx) ## up 

ndata.all <- read.csv("../Data/modTTESTcorALLpvalues_alldatafiltered.txt", sep="\t", comment.char="#", header=TRUE) # 33788    80
annot <- ndata.all[ , c("Symbol", "Entrez_Gene_ID")]

entrzup <- unique(annot[annot$Symbol %in% ssup, ]$Entrez_Gene_ID) ## 58
entrzup <- entrzup[entrzup!=""]
entrzup <- entrzup[!is.na(entrzup)] 
length(entrzup)##57
qup <- rep(1, length(entrzup))
names(qup) <- unlist(lapply(entrzup,function(x) paste("geneid.",x, sep='')))

entrzdn <- unique(annot[annot$Symbol %in% ssdn, ]$Entrez_Gene_ID) ##60
entrzdn <- entrzdn[entrzdn!=""]
entrzdn <- entrzdn[!is.na(entrzdn)] ## 60
length(entrzdn)
qdn <- rep(-1, length(entrzdn))
names(qdn) <- unlist(lapply(entrzdn,function(x) paste("geneid.",x, sep='')))

qq <- c(qup,qdn)
save(qq, file=paste("../Output/q-", sourceName, "-feb9.RData", sep=""))
}
