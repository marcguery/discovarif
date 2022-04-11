#!/bin/bash
#Config for pipeline (Wed 30 Sep 11:12:18 CEST 2020)

configversion="0.0.3"
configrealversion="0.0.3"

##############################PARAMETERS##############################
#Fill the TOBEFILLED parts from the rawdata files
#Customize the paths to your liking as
# long as all rawdata and output are in the same directory

##RAWDATA##
#Path for all rawdata
export DATADIR=TOBEFILLED #Absolute path to the input folder
export REMOTEDATADIR= #Absolute path to a remote data directory. 
                                # Requires a SSH key between the local and remote machines.
                                # Fill this variable if your files are located on a server.
export REMOTEADDRESS= # Address of the remote machine hosting the files.
                        # Fill this variable if your files are located on a server.
#Fasta file
export GENOME="$DATADIR"/TOBEFILLED #Absolute path to fasta reference genome file
export PLOIDY=1 #Ploidy of the genome
export INFOGENOME="$DATADIR"/TOBEFILLED #Absolute path to tabulated file with core annotations of the genome 
                                               #(see Miles 2016 for Plasmodium: 10.1101/gr.203711.115)
export GFF="$DATADIR"/TOBEFILLED #Absolute path to gff features file

export SAMPLEFILE="$DATADIR"/TOBEFILLED #Absolute path to sample tsv file
#FASTQ Read files location
export READDIR="$DATADIR"/TOBEFILLED
#Adapters that should be clipped (Trimmomatic)
export CLIPS="$DATADIR"/TOBEFILLED

##OUTPUT##
#Path for all output
export OUTDIR=TOBEFILLED #Absolute path to the output folder

export REMOTEOUTDIR= #Absolute path to a remote output directory. 
                                # Requires a SSH key between the local and remote machines
                                # Fill this variable if your files are located on a server.

#Quality files
export QUALDIR="$OUTDIR"/fastqc
#Trimmed read files
export TRIMDIR="$OUTDIR"/trim
export TRIMEXT=".paired.gz"
#BAM files and their indexes
export BAMBAIDIR="$OUTDIR"/bambai
export BAMEXT=".dd.sorted.bam"
#Variant files
export VARIANTDIR="$OUTDIR"/variants
#SNP/small INDEL directory
export SNPDIR="$VARIANTDIR"/SNPs-sINDELs
export GVCFDIR="$SNPDIR"/gvcf
#CNV directory
export CNVDIR="$VARIANTDIR"/CNVs
#DELLY directory
export DELLYDIR="$VARIANTDIR"/Others

##############################----------##############################

##############################UTILS##############################
# Provide the path to all executables

export FASTQC=fastqc
export TRIMMOMATIC=TrimmomaticPE
export BWA=bwa
export SAMTOOLS=samtools
export BCFTOOLS=bcftools
export VARIF=varif
export BEDTOOLS=bedtools
export DELLY=delly
export BEDGRAPHTOBIGWIG=bedGraphToBigWig
export PICARD=PicardCommandLine
export GATK=gatk

##############################-----##############################
