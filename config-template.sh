#!/bin/bash
#Config for pipeline (Wed 30 Sep 11:12:18 CEST 2020)

##############################PARAMETERS##############################
#Fill the TOBEFILLED parts

##RAWDATA##
#Path for all rawdata
export DATADIR=TOBEFILLED #Absolute path to the input folder
REMOTEDATADIR= #Absolute path to a remote data directory. 
                                # Requires a SSH key between the local and remote machines.
                                # Fill this variable if your files are located on a server.
REMOTEADDRESS= # Address of the remote machine hosting the files.
                        # Fill this variable if your files are located on a server.
#Fasta file
export GENOME="$DATADIR"/Genome/TOBEFILLED #Absolute path to fasta reference genome file
export PLOIDY=1 #Ploidy of the genome
export INFOGENOME="$DATADIR"/Genome/TOBEFILLED #Absolute path to tabulated file with core annotations of the genome 
                                               #(see Miles 2016 for Plasmodium: 10.1101/gr.203711.115)
export GFF="$DATADIR"/Genome/TOBEFILLED #Absolute path to gff features file
#FASTQ Read files
export READS="$DATADIR"/reads/TOBEFILLED #Wildcard to select all read files (e. g., "$DATADIR"/reads/*.fastq.gz)
#Regular expression to differentiate samples
export SAMPEXP='TOBEFILLED' #Contains all the sequence of characters in read file names identifying samples
                            #(e. g., 'sample[0-9]{1,2}' if the samples are identified as: sample1, sample2, ..., 
                            # sample99')
number_of_samples=TOBEFILLED
#Reads
export IDR1="TOBEFILLED" #Unique sequence of characters in read file names identifying reads 1 (e. g. "R1", or "_R1_")
export IDR2="TOBEFILLED" #Unique sequence of characters in read file names identifying reads 2 (e. g. "R2", or "_R2_")
#Control sample
export CONTROLNAME="TOBEFILLED" #Name of the control sample used for the CNV discovery. It is the name of the
                                #corresponding bam file, without the extensions (e. g., if the bam file is 
                                #"control.sorted.bam", its name is "control")

##OUTPUT##
#Path for all output
export OUTDIR=TOBEFILLED #Absolute path to the output folder

REMOTEOUTDIR= #Absolute path to a remote output directory. 
                                # Requires a SSH key between the local and remote machines
                                # Fill this variable if your files are located on a server.

#Quality files
export QUALDIR="$OUTDIR"/fastqc
export TRIMDIR="$OUTDIR"/trim
#Trimmed files
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

##############################SERVER##############################
#This script will copy your files from a distant server to your current session
# if they are missing from your current session and if you provided a remote address

if [ ! -z $REMOTEADDRESS -a $dry -eq 1 ];then
    echo "Copying missing files in $DATADIR from $REMOTEADDRESS:$REMOTEDATADIR"
    mkdir -p "$DATADIR"
    ssh "$REMOTEADDRESS" [ -d "$REMOTEDATADIR" ] || \
        { echo "$REMOTEDATADIR does not exist in $REMOTEADDRESS"; exit 1; }
    rsync -a --ignore-existing --progress "$REMOTEADDRESS":"$REMOTEDATADIR"/ "$DATADIR"/

    echo "Copying missing files in $OUTDIR from $REMOTEADDRESS:$REMOTEOUTDIR"
    mkdir -p "$OUTDIR"
    ssh "$REMOTEADDRESS" [ -d "$REMOTEOUTDIR" ] || \
        { echo "$REMOTEOUTDIR does not exist in $REMOTEADDRESS"; exit 1; }
    rsync -a --ignore-existing --progress "$REMOTEADDRESS":"$REMOTEOUTDIR"/ "$OUTDIR"/
fi

##############################----------##############################

##############################STATS##############################

SAMPLES=($(ls -1 $READS | sort | grep -o -E $SAMPEXP | uniq))
#Number of samples
export SAMPLENUM=${#SAMPLES[@]}

declare -i oktocontinue=1
[ ! $SAMPLENUM -eq $number_of_samples ] &&
{ echo "Found $SAMPLENUM samples but expected $number_of_samples"; \
echo "The sorted samples"; \
echo "$(ls -1 $READS | sort)"; \
echo "were identified with: ${SAMPLES[@]}"; oktocontinue=1; } ||
{ echo "Successfully indexed $SAMPLENUM samples"; oktocontinue=0; }

[ ! $oktocontinue -eq 0 ] && exit 1

##############################----------##############################

##############################UTILS##############################
# Provide the path to all executables and config files

export FASTQC=fastqc
export TRIMMOMATIC=TrimmomaticPE
export CLIPS="$DATADIR"/adapters/TOBEFILLED
export BWA=bwa
export SAMTOOLS=samtools
export BCFTOOLS=bcftools
export VARIF=varif
export BEDTOOLS=bedtools
export DELLY=delly
export DELLYSAMPLES="$DATADIR"/dellysamples/TOBEFILLED #Tablulated file containing 'tumor' and 'control'
                                                       #samples (see DELLY documentation in somatic mode)
export BEDGRAPHTOBIGWIG=bedGraphToBigWig
export PICARD=PicardCommandLine
export GATK=gatk

##############################-----##############################
