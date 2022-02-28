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
export GENOME="$DATADIR"/TOBEFILLED #Absolute path to fasta reference genome file
export PLOIDY=1 #Ploidy of the genome
export INFOGENOME="$DATADIR"/TOBEFILLED #Absolute path to tabulated file with core annotations of the genome 
                                               #(see Miles 2016 for Plasmodium: 10.1101/gr.203711.115)
export GFF="$DATADIR"/TOBEFILLED #Absolute path to gff features file

export SAMPLEFILE="$DATADIR"/TOBEFILLED #Absolute path to sample tsv file
#FASTQ Read files location
export READDIR="$DATADIR"/TOBEFILLED
#Adapters that shouuld be clipped (Trimmomatic)
export CLIPS="$DATADIR"/TOBEFILLED

##OUTPUT##
#Path for all output
export OUTDIR=TOBEFILLED #Absolute path to the output folder

REMOTEOUTDIR= #Absolute path to a remote output directory. 
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

##############################CHECK##############################

[ ! -f $SAMPLEFILE ] && { "Sample file $SAMPLEFILE does not exist"; exit 1; }
SAMPLES=($(cut -f1 $SAMPLEFILE | tail -n+2))
#Number of samples
export SAMPLENUM=$(($(cut -f1 $SAMPLEFILE | tail -n+2 | sort | uniq | wc -l)))

[ ! $SAMPLENUM -eq ${#SAMPLES[@]} ] &&
{ echo "Found $SAMPLENUM uniquely identified samples but expected ${#SAMPLES[@]}"; \
echo "The sorted samples with their number of occurences (should all be unique):"; \
echo "$(cut -f1 $SAMPLEFILE | tail -n+2 | sort | uniq -c)"; exit 1; }

##############################----------##############################

##############################UTILS##############################
# Provide the path to all executables and config files

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
