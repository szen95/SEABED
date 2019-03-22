library(gplots)

# read file containing cell line information
cell.info <- read.table("./example_inputs/cell_line_info.csv",sep=",", header=T)
cell.info <- cell.info[-which(is.na(cell.info[,"Tissue"])),]

# read a binary mutation event matrix (BEM) for genes X cell lines
mut <- read.table("./example_inputs/PANCAN_simple_MOBEM.tsv", row.names=1, sep="\t", check.names=F, header=T)
mut <- mut[-grep("_HypMET", rownames(mut)),]

# directory containing output files from segmentation algorithm
directory <- "./example_inputs/Segmentation_output"

# add tissue labels to BEM
cell.tissue <- as.character(cell.info$Tissue[match(colnames(mut), cell.info$sample_name)])
features <- c()
for(i in na.omit(unique(cell.tissue))){
  feature.bin <- rep(0, ncol(mut))
  feature.bin[which(cell.tissue == i)] <- 1
  features <- rbind(features, feature.bin)
}
rownames(features) <- na.omit(unique(cell.tissue))
colnames(features) <- colnames(mut)
mut <- rbind(mut, features)

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


subpop.enrich <- function(dat, feature){
  subpops <- unique(dat[,"Subpop."])
  result <- c()
  for(i in subpops){
    yy <- length(which(dat[,feature] == 1 & dat[,"Subpop."] == i))
    yn <- length(which(dat[,feature] == 1 & dat[,"Subpop."] != i))
    ny <- length(which(dat[,feature] != 1 & dat[,"Subpop."] == i))
    nn <- length(which(dat[,feature] != 1 & dat[,"Subpop."] != i))
    p.val <- fisher.test(matrix(c(yy, yn,ny,nn), nrow=2, ncol=2), alternative="greater")$p.value
    prop <- yy/length(which(dat[,"Subpop."] == i)) * 100
    result <- rbind(result, c(i, feature, round(prop,1), p.val))
  }
  colnames(result) <- c("Subpop.", "Feature", "Percent in Subpop.", "P-value")
  return(result)
}


fs <- dir(directory)
fs <- fs[grep("SubpopDataALL.csv", fs)]

for(file.ind in 1:length(fs)){
fc <- file(paste(directory,fs[file.ind], sep=""))
drug.pair <- paste(unlist(strsplit( fs[file.ind], split="_"))[4], collapse="_")

tmp <- strsplit(readLines(fc), ",")
close(fc)
tmp <- tmp[-1]
# clusts <- lapply(tmp, function(x){x[3:(as.numeric(x[2]) +2)]})
clusts <- lapply(tmp, function(x){x[13:(as.numeric(x[12]) +12)]})
names(clusts) <- unlist(lapply(tmp, function(x){x[11]}))


# Get cell subpop classes
cell.clusts <- c()
for(i in 1:length(clusts)){
  cells <- clusts[[i]]
  tmp <- rep(names(clusts)[i], length(cells))
  names(tmp) <- cells
  cell.clusts <- c(cell.clusts, tmp)
}

dat <- matrix(cell.clusts, nrow=length(cell.clusts), ncol=1, dimnames=list(names(cell.clusts), "Subpop."))

# binary encoding for mutation data
samps <- match(rownames(dat), remove.sep(colnames(mut)))
mut.tmp <- mut[,na.omit(samps)]
dat <- dat[which(!is.na(samps)),, drop=F]

dat <- cbind(dat, t(mut.tmp))
class(dat) <- "numeric"

# remove low frequency mutations (<5%)
dat <- dat[,-c(which(apply(dat[,-1],2, sum) < nrow(dat)/20) + 1)]

enrich.tests <- c()
for(feature in colnames(dat[,-1])){
  enrich.tests <- rbind(enrich.tests, subpop.enrich(dat, feature))
}

enrich.tests <- cbind(enrich.tests, p.adjust(as.numeric(enrich.tests[,"P-value"]), method="BH"))
colnames(enrich.tests)[5] <- "Adjusted P-value"
enrich.tests <- enrich.tests[order(as.numeric(enrich.tests[,4])),]

# get numbers of cell lines
cell.count <- sapply(enrich.tests[,1], function(x){length(clusts[[as.character(x)]])})
if(length(cell.count) == 0){
cell.count <- 0
}

feature.count <- round(cell.count * as.numeric(enrich.tests[,3])/100,0)
       
enrich.tests <- cbind(enrich.tests[,1, drop=F], "Number of cells in Subpop." = cell.count, enrich.tests[,2:3], "Number of cells with Feature" = feature.count, enrich.tests[,4:5,drop=F])

# write outputs of biomarkers for each drug pair comparison
write.table(enrich.tests,paste("./Biomarkers/", drug.pair,"_",  ".txt", sep=""), sep="\t",row.names=F, quote=F)
}

