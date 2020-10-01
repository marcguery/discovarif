#!/bin/bash
#Config for pipeline (Wed 30 Sep 11:12:18 CEST 2020)

##############################PARAMETERS##############################
#Fill the TOBEFILLED parts

##RAWDATA##
#Path for all rawdata
export DATADIR=TOBEFILLED
#Fasta file
export GENOME="$DATADIR"/Genome/TOBEFILLED
export INFOGENOME="$DATADIR"/Genome/TOBEFILLED
export GFF="$DATADIR"/Genome/TOBEFILLED
#FASTQ Read files
export READS="$DATADIR"/reads/TOBEFILLED
export SAMPLENUM=$(bc <<< $(($(ls -1 $READS | wc -l)/2)))
echo "Found $SAMPLENUM samples"
#Regular expression to differentiate samples
export SAMPEXP='TOBEFILLED'
SAMPLES=($(ls -1 $READS | grep -o -E $SAMPEXP | uniq))
export IDR1="TOBEFILLED"
export IDR2="TOBEFILLED"
export CONTROLNAME="TOBEFILLED"

##OUTPUT##
#Path for all output
export OUTDIR=TOBEFILLED
#Quality files
export QUALDIR="$OUTDIR"/fastqc
export TRIMDIR="$OUTDIR"/trim
#Trimmed failes
export TRIMREADS="$TRIMDIR"/*.paired.gz
#SAM files
export SAMDIR="$OUTDIR"/sam
#BAM files and their indexes
export BAMBAIDIR="$OUTDIR"/bambai

export BAMFILES="$BAMBAIDIR"/*.sorted.bam
#Variant files
export VARIANTDIR="$OUTDIR"/variants
#SNP/small INDEL directory
export SNPDIR="$VARIANTDIR"/SNPs-sINDELs
#CNV directory
export CNVDIR="$VARIANTDIR"/CNVs
#DELLY directory
export DELLYDIR="$VARIANTDIR"/Others

##############################----------##############################

##############################UTILS##############################

export FASTQC=TOBEFILLED
export TRIMMOMATIC=TOBEFILLED
export CLIPS="$DATADIR"/adapters/TOBEFILLED
export BWA=TOBEFILLED
export SAMTOOLS=TOBEFILLED
export BCFTOOLS=TOBEFILLED
export VARIF=TOBEFILLED
export BEDTOOLS=TOBEFILLED
export DELLY=TOBEFILLED
export DELLYSAMPLES="$DATADIR"/dellysamples/TOBEFILLED

##############################-----##############################
