#!/bin/bash
#Full pipeline (Wed 30 Sep 11:12:18 CEST 2020)
source config.sh

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i doQmv=0
declare -i doSnp=0
declare -i doCnv=0
declare -i doOth=0
declare -i maxthreads=1
while getopts ":mscot:h" o; do
    case "${o}" in
        m) # Launch quality/mapping/mpileup step (NOT AVAILABLE).
            doQmv=1
            usage
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
        t) # Launch this number of processes in parallel
            maxthreads=$((OPTARG))
            ;;
        h | *) # Show help.
            usage
            ;;
    esac
done
shift $((OPTIND-1))
[ $maxthreads -gt $SAMPLENUM ] && { maxthreads=$SAMPLENUM; }
##############################-------##############################

##############################PIPELINE##############################

####Quality and alignment with each sample####
if [ $doQmv -eq 1 ];then
    echo "Doing the Quality/mapping/variant step for each sample..."
    numsamples=$(bc <<< $SAMPLENUM/$maxthreads)
    number=0
    while [ $number -lt $(($SAMPLENUM-$numsamples)) ];do
        nextid=$(($number+$numsamples))
        echo "New batch: ${SAMPLES[$number]} to ${SAMPLES[$((nextid-1))]}"
        ./mapper-caller.sh -qm -s $number:$nextid &
        number=$nextid
    done
    echo "New batch: ${SAMPLES[$number]} to ${SAMPLES[$(($SAMPLENUM-1))]}"
    ./mapper-caller.sh -qm -s $number:$SAMPLENUM &
    wait
    echo "Done the Quality/mapping/variant step for each sample!"
fi
########

####SNP/small INDEL filtering####
if [ $doSnp -eq 1 ];then
    echo "Doing the SNPs/INDELs step..."
    mkdir -p "$SNPDIR"

    #Doing the mpileup with all the samples at once
    $BCFTOOLS mpileup -f $GENOME $BAMFILES -Ou -d 99999 -a DP,AD,SP | \
    $BCFTOOLS call -m -v -Ob -o "$SNPDIR"/all-samples.bcf

    $BCFTOOLS view "$SNPDIR"/all-samples.bcf > "$SNPDIR"/all-samples.vcf
    $BCFTOOLS view -i '%QUAL>=10' "$SNPDIR"/all-samples.bcf > "$SNPDIR"/all-samples-Qual10.vcf

    ##Filtering variants
    #In renamed vcf files, sample names replace file paths
    $VARIF -vcf "$SNPDIR"/all-samples-renamed.vcf -gff "$GFF" -fasta "$GENOME" \
    --no-fixed --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.6 \
    --ratio-no-alt 0.4 --csv "$SNPDIR"/filtered-SNPs-sINDELs.csv \
    --filteredvcf "$SNPDIR"/filtered-SNPs-sINDELs.vcf

    $VARIF -vcf "$SNPDIR"/all-samples-Qual10-renamed.vcf -gff "$GFF" -fasta "$GENOME" \
    --no-fixed --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.6 \
    --ratio-no-alt 0.4 --csv "$SNPDIR"/filtered-SNPs-sINDELs-Qual10.csv \
    --filteredvcf "$SNPDIR"/filtered-SNPs-sINDELs-Qual10.vcf
    echo "Done the SNPs/INDELs step!"
fi
########

####CNV####
if [ $doCnv -eq 1 ];then
    # ##Setup
    # echo "Doing the CNVs step..."
    # mkdir -p "$CNVDIR"/view

    # echo "###Background preparation###"
    # echo "Getting chromosome sizes..."
    # $SAMTOOLS faidx "$GENOME"
    # cut -f1,2 "$GENOME".fai > "$CNVDIR"/chrom.sizes
    # echo "Getting core genome..."
    # grep "Core" "$INFOGENOME" | cut -f1-3 > "$CNVDIR"/view/3D7-core.bed

    # echo "Getting CDS coordinates..."
    # grep CDS $GFF | cut -f1,4,5 | sort -V > "$CNVDIR"/view/cds.bed
    # $BEDTOOLS intersect -a "$CNVDIR"/view/3D7-core.bed -b "$CNVDIR"/view/cds.bed > "$CNVDIR"/view/cds-core.bed

    # echo "Getting CDS perbase location..."
    # $BEDTOOLS genomecov -g "$CNVDIR"/chrom.sizes -i "$CNVDIR"/view/cds.bed -dz | cut -f1,2 > "$CNVDIR"/view/cds-perbase.bed
    
    # ##Obtaining cov files
    # numsamples=$(bc <<< $SAMPLENUM/$maxthreads)
    # number=0
    # while [ $number -lt $(($SAMPLENUM-$numsamples)) ];do
    #     nextid=$(($number+$numsamples))
    #     echo "New batch: ${SAMPLES[$number]} to ${SAMPLES[$((nextid-1))]}"
    #     ./get-coverage.sh -s $number:$nextid &
    #     number=$nextid
    # done
    # echo "New batch: ${SAMPLES[$number]} to ${SAMPLES[$(($SAMPLENUM-1))]}"
    # ./get-coverage.sh -s $number:$SAMPLENUM &
    # wait
    
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
