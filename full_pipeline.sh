#!/bin/bash
#Full pipeline (Wed 30 Sep 11:12:18 CEST 2020)
source config.sh

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i doQma=0
declare -i doSnp=0
declare -i doCnv=0
declare -i doOth=0
declare -i maxthreads=1
declare -i maxthreadspersample=1
while getopts ":mscot:u:h" o; do
    case "${o}" in
        m) # Launch quality/mapping/mpileup step.
            doQma=1
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
        u) # Launch this number of processes in parallel
            maxthreadspersample=$((OPTARG))
            ;;
        h | *) # Show help.
            usage
            ;;
    esac
done
shift $((OPTIND-1))
[ $maxthreadspersample -gt $maxthreads ] && \
{ echo "Error: Maximum allowed threads of $maxthreads while $maxthreadspersample requested"; exit 1; }
threadsamples=$(bc <<< $maxthreads/$maxthreadspersample)
[ $threadsamples -gt $SAMPLENUM ] && { threadsamples=$SAMPLENUM; }
##############################-------##############################

##############################PIPELINE##############################

####Quality and alignment with each sample####
if [ $doQma -eq 1 ];then
    echo "Doing the Quality/mapping/variant step for each sample..."
    start_index=0
    number_of_threads=7
    sequential_per_process=0
    iter=0
    echo ${SAMPLES[@]}
    while [ $start_index -le $(($SAMPLENUM-$sequential_per_process)) ];do
        sequential_per_process=$(bc <<< $(($SAMPLENUM-$start_index))/$(($threadsamples-$iter)))
        end_index=$(($start_index+$sequential_per_process))
        echo "New batch: ${SAMPLES[$start_index]} to ${SAMPLES[$(($end_index-1))]}"
        echo "./mapper-caller.sh -m -s $start_index:$end_index -t $maxthreadspersample"
        start_index=$end_index
        ((iter++))
    done
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
    $BCFTOOLS call --threads $(($maxthreads-1)) -m -v -Ob -o "$SNPDIR"/all-samples.bcf

    $BCFTOOLS view "$SNPDIR"/all-samples.bcf > "$SNPDIR"/all-samples.vcf
    $BCFTOOLS view -i '%QUAL>=10' "$SNPDIR"/all-samples.bcf > "$SNPDIR"/all-samples-Qual10.vcf

    ##Filtering variants
    #In renamed vcf files, sample names replace file paths
    $VARIF -vcf "$SNPDIR"/all-samples-renamed.vcf -gff "$GFF" -fasta "$GENOME" \
    --no-fixed --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.6 \
    --ratio-no-alt 0.4 --csv "$SNPDIR"/alt06ref04/filtered-SNPs-sINDELs-0604.csv \
    --filteredvcf "$SNPDIR"/alt06ref04/filtered-SNPs-sINDELs-0604.vcf

    $VARIF -vcf "$SNPDIR"/all-samples-Qual10-renamed.vcf -gff "$GFF" -fasta "$GENOME" \
    --no-fixed --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.6 \
    --ratio-no-alt 0.4 --csv "$SNPDIR"/alt06ref04/filtered-SNPs-sINDELs-Qual10-0604.csv \
    --filteredvcf "$SNPDIR"/alt06ref04/filtered-SNPs-sINDELs-Qual10-0604.vcf

    $VARIF -vcf "$SNPDIR"/all-samples-renamed.vcf -gff "$GFF" -fasta "$GENOME" \
    --no-fixed --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.55555 \
    --ratio-no-alt 0.05 --csv "$SNPDIR"/alt05ref005/filtered-SNPs-sINDELs-05005.csv \
    --filteredvcf "$SNPDIR"/alt06ref005/filtered-SNPs-sINDELs-05005.vcf

    $VARIF -vcf "$SNPDIR"/all-samples-renamed.vcf -gff "$GFF" -fasta "$GENOME" \
    --all-variants --all-regions --no-show --depth 6 --ratio-alt 0.8 \
    --ratio-no-alt 0.02 --csv "$SNPDIR"/alt08ref002/filtered-SNPs-sINDELs-08002-all.csv \
    --filteredvcf "$SNPDIR"/alt08ref002/filtered-SNPs-sINDELs-08002-all.vcf
    echo "Done the SNPs/INDELs step!"
fi
########

####CNV####
if [ $doCnv -eq 1 ];then
    ##Setup
    echo "Doing the CNVs step..."
    mkdir -p "$CNVDIR"/view

    echo "###Background preparation###"
    echo "Getting chromosome sizes..."
    $SAMTOOLS faidx "$GENOME"
    cut -f1,2 "$GENOME".fai > "$CNVDIR"/chrom.sizes
    echo "Getting core genome..."
    grep "Core" "$INFOGENOME" | cut -f1-3 > "$CNVDIR"/view/3D7-core.bed

    echo "Getting CDS coordinates..."
    grep CDS $GFF | cut -f1,4,5 | sort -V > "$CNVDIR"/view/cds.bed
    $BEDTOOLS intersect -a "$CNVDIR"/view/3D7-core.bed -b "$CNVDIR"/view/cds.bed > "$CNVDIR"/view/cds-core.bed

    echo "Getting CDS perbase location..."
    $BEDTOOLS genomecov -g "$CNVDIR"/chrom.sizes -i "$CNVDIR"/view/cds.bed -dz | cut -f1,2 > "$CNVDIR"/view/cds-perbase.bed
    
    ##Obtaining cov files
    start_index=0
    number_of_threads=7
    sequential_per_process=0
    iter=0
    echo ${SAMPLES[@]}
    while [ $start_index -le $(($SAMPLENUM-$sequential_per_process)) ];do
        sequential_per_process=$(bc <<< $(($SAMPLENUM-$start_index))/$(($threadsamples-$iter)))
        end_index=$(($start_index+$sequential_per_process))
        echo "New batch: ${SAMPLES[$start_index]} to ${SAMPLES[$(($end_index-1))]}"
        ./get-coverage.sh -s $start_index:$end_index &
        start_index=$end_index
        ((iter++))
    done
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
    ./DELLY-caller.sh -s "$(ls -1 $BAMFILES | grep -v -E $CONTROLNAME)" -c "$(ls -1 $BAMFILES | grep -E $CONTROLNAME)"
    echo "Done the Other variants step!"
fi
########

##############################--------##############################
