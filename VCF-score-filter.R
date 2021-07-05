#!/usr/bin/env Rscript
#Filter VCF entries based on their score variants (Mon  5 Jul 11:09:22 CEST 2021)
if (as.integer(R.Version()$major) < 3) stop("Wrong version number, R 3 at least")

##############################LIBRARIES##############################

##############################---------##############################

##############################PARAMETERS##############################

args <- commandArgs(trailingOnly=T)
options(stringsAsFactors = F)

vcffilename <- args[1] #Name of the input VCF file
propbadentries <- args[2] #Proportion of estimated false positive
outfilename <- args[3] #Name of the filtered output VCF file

##############################----------##############################

##############################PIPELINE##############################

##### Read twice the vcf file, first for the columns names, second for the data
whole_vcf<-readLines(vcffilename)
vcf_data<-read.table(vcffilename, stringsAsFactors = FALSE)

#### Filter for the columns names
vcf <- whole_vcf[-(grep("#CHROM",whole_vcf)+1):-(length(whole_vcf))]
vcf_names <- unlist(strsplit(vcf[length(vcf)],"\t"))
names(vcf_data) <- vcf_names

#### Calculating the score cutoff
scores <- as.numeric(vcf_data$QUAL)
cutoff <- quantile(scores, as.numeric(propbadentries))

#### Writing the filtered VCF file
vcfheader <- whole_vcf[grepl("#", whole_vcf)]
write(vcfheader,
      outfilename, 
      sep = "\t")
write.table(vcf_data[vcf_data$QUAL>cutoff,], 
            outfilename,
            append = TRUE, row.names = F, col.names = F, quote = F, sep = "\t")

##############################--------##############################
