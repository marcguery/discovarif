#!/bin/bash
#Full pipeline (Wed 30 Sep 11:12:18 CEST 2020)

source config.sh

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i doQmv=0
declare -i doSnp=0
declare -i doCnv=0
declare -i doOth=0
while getopts ":mscoh" o; do
    case "${o}" in
        m) # Launch quality/mapping/mpileup step.
            doQmv=1
            ;;
        s) # Launch SNP/small INDEL step.
            doSnp=1
            ;;
        c) # Launch CNV step.
            doCnv=1
            ;;
        o) # Launch other variants step.
            doOth=1
            ;;
        h | *) # Show help.
            usage
            ;;
    esac
done
shift $((OPTIND-1))
##############################-------##############################

##############################PIPELINE##############################

####Quality and alignment with each sample####
if [ $doQmv -eq 1 ];then
    echo "Doing the Quality/mapping/variant step for each sample..."
    ./mapper-caller.sh -qm -s 0:2 &
    ./mapper-caller.sh -qm -s 2:4 &
    ./mapper-caller.sh -qm -s 4:6 &
    ./mapper-caller.sh -qm -s 6:8 &
    ./mapper-caller.sh -qm -s 8:10 &
    ./mapper-caller.sh -qm -s 10:12 &
    ./mapper-caller.sh -qm -s 12:14 &
    wait
    echo "Done the Quality/mapping/variant step for each sample!"
fi
########

####SNP/small INDEL filtering####
if [ $doSnp -eq 1 ];then
    echo "Doing the SNPs/INDELs step..."
    mkdir -p "$SNPDIR"

    #Doing the mpileup with all the samples at once
    $BCFTOOLS mpileup -f $GENOME "$BAMBAIDIR"/*.sorted.bam -Ou -d 99999 -a DP,AD,SP | \
    $BCFTOOLS call -m -v -Ob -o "$SNPDIR"/all-samples.bcf

    $BCFTOOLS view -i '%QUAL>=10' "$SNPDIR"/all-samples.bcf > "$SNPDIR"/all-samples-Qual10.vcf

    ##Filtering variants
    $VARIF -vcf "$SNPDIR"/all-samples-Qual10.vcf -gff $GFF -fasta $GENOME \
    --no-fixed --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.8 \
    --ratio-no-alt 0.2 --csv "$SNPDIR"/filtered-SNPs-sINDELs.csv
    echo "Done the SNPs/INDELs step!"
fi
########

####CNV####
if [ $doCnv -eq 1 ];then
    ##Setup
    echo "Doing the CNVs step..."
    mkdir -p "$CNVDIR"/bed

    echo "###Background preparation###"
    echo "Getting chromosome sizes..."
    $SAMTOOLS faidx $GENOME
    cut -f1,2 $GENOME.fai > "$CNVDIR"/chrom.sizes
    echo "Getting core genome..."
    grep "Core" $INFOGENOME | cut -f1-3 > "$CNVDIR"/bed/3D7-core.bed

    echo "Getting CDS coordinates..."
    grep CDS $GFF | cut -f1,4,5 | sort -V > "$CNVDIR"/bed/cds.bed
    $BEDTOOLS intersect -a "$CNVDIR"/bed/3D7-core.bed -b "$CNVDIR"/bed/cds.bed > "$CNVDIR"/bed/cds-core.bed

    echo "Getting CDS perbase location..."
    $BEDTOOLS genomecov -g "$CNVDIR"/chrom.sizes -i "$CNVDIR"/bed/cds.bed -dz | cut -f1,2 > "$CNVDIR"/bed/cds-perbase.bed
    
    ##Obtaining cov files
    ./get-cov.sh -s 0:2 &
    ./get-cov.sh -s 2:4 &
    ./get-cov.sh -s 4:6 &
    ./get-cov.sh -s 6:8 &
    ./get-cov.sh -s 8:10 &
    ./get-cov.sh -s 10:12 &
    ./get-cov.sh -s 12:14 &
    wait
    
    ##Getting filtered CNVs
    mkdir -p "$CNVDIR"/summary
    ./CNV-caller.R "$CNVDIR" "-perbasecds-core.coverage" "$CNVDIR"/summary "-core-cov.tsv" "CNV_withoutAPI-MT.csv" $CONTROLNAME
    
    echo "Done the CNVs step!"
fi
########

####Other variants####
if [ $doOth -eq 1 ];then
    echo "Doing the Other variants step..."
    mkdir -p "$DELLYDIR"
    ./DELLY-caller.sh
    echo "Done the Other variants step!"
fi
########

##############################--------##############################
