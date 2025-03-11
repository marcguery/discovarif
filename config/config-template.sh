#!/bin/bash
#Config for pipeline (Wed 30 Sep 11:12:18 CEST 2020)

configcompversion="0.0.7"
configversion="0.0.6"

##############################PARAMETERS##############################
#Fill the TOBEFILLED parts 
# and optionally the IFAVAILABLE parts

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
export INFOGENOME="$DATADIR"/IFAVAILABLE #Absolute path to tabulated file with core annotations of the genome 
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
#Intermediate mapping files location (will be deleted at the end of the run)
export TMPSAMDIR="$BAMBAIDIR"/_tmp
#Variant files
export VARIANTDIR="$OUTDIR"/variants
#SNP/small INDEL directory
export SNPDIR="$VARIANTDIR"/SNPs-sINDELs
#Intermediate GVCF files location (will be deleted at the end of the run)
export TMPGVCFDIR="$SNPDIR"/_tmp
#GATK's GVCFs that were successfully created
export GVCFDIR="$SNPDIR"/gvcf
#CNV directory
export CNVDIR="$VARIANTDIR"/CNVs
#DELLY directory
export DELLYDIR="$VARIANTDIR"/Others

##############################----------##############################

##############################UTILS##############################
# Provide the path to all executables

export FASTQC=fastqc
export BBDUK=bbduk.sh
export BWA=bwa
export SAMTOOLS=samtools
export BCFTOOLS=bcftools
export ALFRED=alfred
export VARIF=varif
export BEDTOOLS=bedtools
export DELLY=delly
export BEDGRAPHTOBIGWIG=bedGraphToBigWig
export GATK=gatk
export SEQKIT=seqkit

# Provide the path to all JAR files
export PICARD_jar=/usr/share/java/picard.jar

##############################-----##############################
