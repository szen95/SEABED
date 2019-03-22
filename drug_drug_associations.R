# Functions to classify drug pairs based on the distribution of subpopulations
library(gplots)
setwd("~/Google Drive/AZ/POET/POET-manuscript/Methods/")

# remove separaters in cell names
remove.sep <- function(query.name){
  query.name <- gsub("_", "", query.name)
  query.name <- gsub("\\.", "", query.name)
  query.name <- gsub("\\(", "", query.name)
  query.name <- gsub("\\)", "", query.name)
  query.name <- gsub("-", "", query.name)
  query.name <- gsub(" ", "", query.name)
  return(toupper(query.name))
}

subpop.position <- function(subpop.dat, subpop){
  dA = subpop.dat$`Drug A Average IC50`[subpop]
  dB = subpop.dat$`Drug B Average IC50`[subpop]
  thresA = subpop.dat$`Drug A IC50 20% cutoff`[1]
  thresB = subpop.dat$`Drug B IC50 20% cutoff`[1]
  rangeA = abs(subpop.dat$`Drug A max IC50 conc.`[1] - subpop.dat$`Drug A min IC50 conc.`[1])
  rangeB = abs(subpop.dat$`Drug B max IC50 conc.`[1] - subpop.dat$`Drug B min IC50 conc.`[1])
  result = rep(NA, length(dA))
  result[which(dA < thresA & dB < thresB)] = "both sensitive"
  result[which(dA < thresA & dB >= thresB)] = "drugA sensitive"
  result[which(dA >= thresA & dB < thresB)] = "drugB sensitive"
  result[which(dA >= thresA & dB >= thresB)] = "both resistant"
  thres.distA <- abs(dA - thresA)
  thres.distB <- abs(dB - thresB)
  dist.normA <- abs(dA - thresA) / rangeA
  dist.normB <- abs(dB - thresB) / rangeB
  return(cbind("Subpop"=result, "drugA distance"= thres.distA, "drugB distance"= thres.distB, "drugA distance normalized"= dist.normA, "drugB distance normalized"= dist.normB))
}

directory <- "./example_inputs/Segmentation_output/"
fs <- dir(directory)
fs <- fs[grep("ExportedSubpopDataALL.csv", fs)]

drugA <- rep(NA, length(fs))
drugB <-  rep(NA, length(fs))
category <-  rep(NA, length(fs))
for(file.ind in 1:length(fs)){

  tmp <- read.table(paste(directory,fs[file.ind], sep=""), sep=",", header=T, check.names = F)
  
  drugA[file.ind] <- as.character(tmp[1,"Drug A"])
  drugB[file.ind] <-  as.character(tmp[1, "Drug B"])
  drugA.cutoff <- tmp[1,]$`Drug A IC50 20% cutoff`
  drugB.cutoff <- tmp[1,]$`Drug B IC50 20% cutoff`
  
  mat <- matrix(0, nrow=2, ncol=2, dimnames=list(c("DrugA.sens", "DrugA.res"), c("DrugB.sens","DrugB.res")))

  mat[1,1] <- sum(tmp[which(tmp$`Drug A Average IC50` < drugA.cutoff & tmp$`Drug B Average IC50` < drugB.cutoff ), "Num of cell lines"])
  mat[1,2] <- sum(tmp[which(tmp$`Drug A Average IC50` < drugA.cutoff & tmp$`Drug B Average IC50` >= drugB.cutoff ), "Num of cell lines"])
  mat[2,1] <- sum(tmp[which(tmp$`Drug A Average IC50` >= drugA.cutoff & tmp$`Drug B Average IC50` < drugB.cutoff ), "Num of cell lines"])
  mat[2,2] <- sum(tmp[which(tmp$`Drug A Average IC50` >= drugA.cutoff & tmp$`Drug B Average IC50` >= drugB.cutoff ), "Num of cell lines"])
  

  # preferential for Drug A
  if( mat[1,2] + mat[2,1] != 0 & mat[1,2] + mat[2,2] != 0){
    if( binom.test(mat[1,2], mat[1,2] + mat[2,1], 0.5, alternative="greater")$p.value < 0.05 & binom.test(mat[1,2], mat[1,2] + mat[2,2], 0.2, alternative="greater")$p.value < 0.05){
      category[file.ind] <- 2
    }
  }
  # preferential for Drug B
  if( mat[1,2] + mat[2,1] != 0 & mat[2,2] + mat[2,1] != 0){
    if( binom.test(mat[2,1], mat[1,2] + mat[2,1], 0.5, alternative="greater")$p.value < 0.05 & binom.test(mat[2,1], mat[2,2] + mat[2,1], 0.2, alternative="greater")$p.value < 0.05){
      category[file.ind] <- 3
    }
  }
  # anti-symmetric response
  if( mat[1,2] + mat[2,2] != 0 & mat[2,1] + mat[2,2] !=0 ){
    if( binom.test(mat[1,2], mat[1,2] + mat[2,2], 0.2, alternative="greater")$p.value < 0.05 & binom.test(mat[2,1], mat[2,1] + mat[2,2], 0.2, alternative="greater")$p.value < 0.05){
      category[file.ind] <- 4
    }
  }
  # sensitive to both Drug A and Drug B
  if( mat[2,2] + mat[1,1] != 0 & mat[1,2] + mat[2,1] + mat[1,1] !=0 ){
    if(binom.test(mat[1,1], mat[2,2] + mat[1,1], 0.04, alternative="greater")$p.value < 0.05 & binom.test(mat[1,1], mat[1,2] + mat[2,1] + mat[1,1], 0.25, alternative="greater")$p.value < 0.05){

      category[file.ind] <- 1
    }
  }
  # no preference
  if(is.na(category[file.ind])){
    category[file.ind] <- 0
  }
}

# Make a matrix of drugA X drugB with their relationship categories  
cat.mat <- matrix(NA, nrow=length(unique(drugA)), ncol=length(unique(drugB)), dimnames=list(unique(drugA), unique(drugB)))
for(i in 1:length(drugA)){
  cat.mat[drugA[i], drugB[i]] <- category[i]
}


# Annotate matrix with their gene targets
drug.info <- read.table("../Data/drugInfo.csv", sep=",", header=T)
drugA.targets <- as.character(drug.info$TARGET[match(rownames(cat.mat), drug.info$DRUG.NAME)])
drugB.targets <- as.character(drug.info$TARGET[match(colnames(cat.mat), drug.info$DRUG.NAME)])
cat.mat <- cat.mat[order(drugA.targets),]
cat.mat <- cat.mat[,order(drugB.targets)]
rownames(cat.mat) <- paste(rownames(cat.mat)," (", sort(drugA.targets), ")", sep="")
colnames(cat.mat) <- paste(colnames(cat.mat)," (", sort(drugB.targets), ")", sep="")

# remove selumetinib with afatinib comparisons IF you want to compare MAPK vs AKT pathways only
cat.mat <- cat.mat[, -grep("selumetinib", colnames(cat.mat))]
cat.mat <- cat.mat[-grep("Afatinib", rownames(cat.mat)),]

tmp<- c()
for(i in 1:nrow(cat.mat))
{
  for(j in 1:ncol(cat.mat)){
    tmp<- rbind(tmp, c(rownames(cat.mat)[i], colnames(cat.mat)[j]))
  }
}

# output the list of drug-drug pairs compared
write.table(tmp, "drug_pairs_list.csv", quote=F, row.names=T, col.names=F, sep=",")

# write matrix describing drug-drug relationships
write.table(cat.mat, "drug-drug_relationships.csv", sep=",", col.names=NA)


