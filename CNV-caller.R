#!/usr/local/R-3.6.0/bin/Rscript
#Filter CNV variants (Thu  1 Oct 10:14:38 CEST 2020)

##############################LIBRARIES##############################

library(reshape2)

##############################---------##############################

##############################PARAMETERS##############################

args <- commandArgs(trailingOnly=T)
options(stringsAsFactors = F)

indir <- args[1]
incoveragepattern <- args[2]
outdir <- args[3]
outcovpattern <- args[4]
outsummary <- args[5]
controlsample <- args[6]

##############################----------##############################

##############################PIPELINE##############################

coveragefiles <- dir(indir, pattern = incoveragepattern)
ratiocds.05 <- data.frame()

for (f in coveragefiles){
    ####Summary per sample####
    covname <- sub(incoveragepattern, outcovpattern, f)
    if (file.exists(paste(outdir, covname, sep="/"))){
        print(paste("Skipping", covname, "because file exists"))
        cov.cds <- read.table(paste(outdir, covname, sep="/"), h=F)[,c(2,3,4,7)]
    }else {
        samplefile <- read.table(paste(indir, f, sep="/"))
        samplefile$V5 <- samplefile$V5/median(samplefile$V5)
        cov.cds <- data.frame(aggregate(samplefile$V5,by=list(samplefile$V1, samplefile$V2, samplefile$V3),
                                    function(x) {round(summary(x), 3)}))
        cov.cds <- cov.cds[with(cov.cds, order(Group.1, Group.2, Group.3)), ]
        write.table(cov.cds, paste(outdir, covname, sep="/"), 
                    quote = F, sep="\t", col.names = F)
    }
    ########

    ####Merging all samples####
    mainname <- sub(incoveragepattern, "", f)
    mainname <- make.names(mainname)
    colnames(cov.cds)[ncol(cov.cds)] <- mainname
    if (all(dim(ratiocds.05)==0)){
        ratiocds.05 <- cov.cds
    }else {
        ratiocds.05 <- merge(ratiocds.05, cov.cds)
    }
    ########

}

#####Best CNVs####
score_from_cov <- function(cdsRow, control){
  if (length(which(cdsRow >= 1.8))>=1 & cdsRow[control] < 1.5){
    multFac <- 1
  }
  else if (length(which(cdsRow <= (0.2)))>=1 & cdsRow[control] > 0.5){
    multFac <- -1
  }
  else{
    score <- 0
    return(score)
  }
  score <- multFac*min(100,max(cdsRow)/min(cdsRow))
  return(score)
}
controlnumber <- which(colnames(ratiocds.05) == make.names(controlsample))
print(paste("The control sample is", controlsample, "(", controlnumber, ")"))
ratiocds.05$score <- apply(ratiocds.05[,4:(length(coveragefiles)+3)], 1,
                           FUN = score_from_cov, control=controlnumber-3)

ratiocds.05 <- ratiocds.05[with(ratiocds.05, order(V2, V3, V4)), ]
ratiocds.05 <- ratiocds.05[ratiocds.05$score!=0,]
write.csv(ratiocds.05, paste(outdir, outsummary, sep="/"), 
            quote = F, row.names = F)
########

##############################--------##############################
