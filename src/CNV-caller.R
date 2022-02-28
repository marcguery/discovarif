#!/usr/bin/env Rscript
#Filter CNV variants (Thu  1 Oct 10:14:38 CEST 2020)
if (as.integer(R.Version()$major) < 3) stop("Wrong version number, R 3 at least")
##############################LIBRARIES##############################

library(reshape2)

##############################---------##############################

##############################PARAMETERS##############################

args <- R.utils::commandArgs(trailingOnly=T, asValues=T)
options(stringsAsFactors = F)
indir <- args[1][[1]]
incoveragepattern <- args[2][[1]]
outdir <- args[3][[1]]
outcovpattern <- args[4][[1]]
outsummary <- args[5][[1]]
controlsamples <- strsplit(args[6][[1]], split = " ")[[1]]
ratiotumor <- as.numeric(args[7][[1]])
ratiocontrol <- as.numeric(args[8][[1]])

##############################----------##############################

##############################PIPELINE##############################

coveragefiles <- dir(indir, pattern = incoveragepattern)
ratiocds.05 <- data.frame()

for (f in coveragefiles){
    ####Summary per sample####
    print(f)
    covname <- sub(incoveragepattern, outcovpattern, f)
    if (file.exists(paste(outdir, covname, sep="/"))){
        print(paste("Skipping", covname, "because file exists"))
        cov.cds <- read.table(paste(outdir, covname, sep="/"), h=F)[,c(2,3,4,7)]
    }else {
        samplefile <- read.table(paste(indir, f, sep="/"))
        samplemedian <- max(median(samplefile$V5), 1)
        if (samplemedian <= 1){
          print(paste("Very low coverage for sample", f, "Results are not reliable"))
        }
        samplefile$V5 <- samplefile$V5/samplemedian
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

}

#####Best CNVs####
score_from_cov <- function(cdsRow, control, mintumor, mincontrol){
  if (length(which(cdsRow >= 1.8)) >= mintumor & length(which(cdsRow[control] < 1.5)) >= mincontrol ){
    multFac <- 1
  }
  else if (length(which(cdsRow <= (0.05))) >= mintumor & length(which(cdsRow[control] > 0.5)) >= mincontrol){
    multFac <- -1
  }
  else{
    score <- 0
    return(score)
  }
  score <- multFac*min(100,max(cdsRow)/min(cdsRow))
  return(score)
}
if (! all(make.names(controlsamples)%in%colnames(ratiocds.05))){
  write.csv(ratiocds.05, paste(outdir, outsummary, sep="/"), 
            quote = F, row.names = F)
  print(colnames(ratiocds.05[,4:(length(coveragefiles)+3)]))
  stop(paste0("The control samples \"", 
  controlsamples, 
  "\" were not all not found in the dataset above, aborting..."))
}

controlindexes <- which(colnames(ratiocds.05) %in% make.names(controlsamples))
mincontrol <- round(ceiling(length(controlindexes)/(1/ratiocontrol)))
tumornumber <- length(colnames(ratiocds.05)) - length(controlindexes) - 3
mintumor <- round(ceiling(tumornumber/(1/ratiotumor)))

print(paste0("The control sample ", as.character(controlsamples), " is at the ", as.character(controlindexes), "th column"))
print(paste0("Looking for variants present in at least ", mintumor," of the ", tumornumber
  ," samples and absent in at least ", mincontrol, " of the ",
  length(controlindexes), " control samples..."))
ratiocds.05$score <- apply(ratiocds.05[,4:(length(coveragefiles)+3)], 1,
                           FUN = score_from_cov, control=controlindexes-3, mintumor=mintumor, mincontrol=mincontrol)

ratiocds.05 <- ratiocds.05[with(ratiocds.05, order(V2, V3, V4)), ]
ratiocds.05 <- ratiocds.05[ratiocds.05$score!=0,]
write.csv(ratiocds.05, paste(outdir, outsummary, sep="/"), 
            quote = F, row.names = F)
########

##############################--------##############################
