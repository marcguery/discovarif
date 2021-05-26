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
        u) # Launch this number of processes for each sample
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
    $BWA index "$GENOME"

    start_index=0
    sequential_per_process=0
    iter=0
    echo ${SAMPLES[@]}
    while [ $start_index -le $(($SAMPLENUM-$sequential_per_process)) ];do
        sequential_per_process=$(bc <<< $(($SAMPLENUM-$start_index))/$(($threadsamples-$iter)))
        end_index=$(($start_index+$sequential_per_process))
        echo "New batch: ${SAMPLES[$start_index]} to ${SAMPLES[$(($end_index-1))]}"
        ./mapper-caller.sh -m -s $start_index:$end_index -t $maxthreadspersample &
        start_index=$end_index
        ((iter++))
    done
    wait
    echo "Done the Mapping step for each sample!"
fi
########

####SNP/small INDEL filtering####
if [ $doSnp -eq 1 ];then
    echo "Doing the SNPs/INDELs step..."
    mkdir -p "$SNPDIR"

    $GATK CreateSequenceDictionary -R $GENOME

    start_index=0
    sequential_per_process=0
    iter=0
    echo ${SAMPLES[@]}
    while [ $start_index -le $(($SAMPLENUM-$sequential_per_process)) ];do
        sequential_per_process=$(bc <<< $(($SAMPLENUM-$start_index))/$(($threadsamples-$iter)))
        end_index=$(($start_index+$sequential_per_process))
        echo "New batch: ${SAMPLES[$start_index]} to ${SAMPLES[$(($end_index-1))]}"
        ./mapper-caller.sh -v -s $start_index:$end_index -t $maxthreadspersample &
        start_index=$end_index
        ((iter++))
    done
    wait
    echo "Done the variant step for each sample!"

    $GATK CombineGVCFs \
        -R $GENOME \
        --variant $(sed -e 's/ / --variant /g' <(echo $SNPDIR/*.vcf.gz)) \
        -O $SNPDIR/cohort.g.vcf.gz

    $GATK --java-options "-Xmx4g" GenotypeGVCFs \
        -R $GENOME \
        -V $SNPDIR/cohort.g.vcf.gz \
        -O $SNPDIR/output.vcf.gz
    
    mkdir -p "$SNPDIR"/alt08ref02
    $VARIF -vcf <(gunzip -c $SNPDIR/output.vcf.gz) -gff "$GFF" -fasta "$GENOME" \
    --no-fixed --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.8 \
    --ratio-no-alt 0.2 --csv "$SNPDIR"/alt08ref02/filtered-SNPs-sINDELs-0802.csv \
    --filteredvcf "$SNPDIR"/alt08ref02/filtered-SNPs-sINDELs-0802.vcf

    ##Filtering variants
    #In renamed vcf files, sample names replace file paths
    mkdir -p "$SNPDIR"/alt06ref04
    $VARIF -vcf <(gunzip -c "$SNPDIR"/output.vcf.gz) -gff "$GFF" -fasta "$GENOME" \
    --no-fixed --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.6 \
    --ratio-no-alt 0.4 --csv "$SNPDIR"/alt06ref04/filtered-SNPs-sINDELs-0604.csv \
    --filteredvcf "$SNPDIR"/alt06ref04/filtered-SNPs-sINDELs-0604.vcf

    mkdir -p "$SNPDIR"/alt05ref005
    $VARIF -vcf <(gunzip -c "$SNPDIR"/output.vcf.gz) -gff "$GFF" -fasta "$GENOME" \
    --no-fixed --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.55555 \
    --ratio-no-alt 0.05 --csv "$SNPDIR"/alt05ref005/filtered-SNPs-sINDELs-05005.csv \
    --filteredvcf "$SNPDIR"/alt05ref005/filtered-SNPs-sINDELs-05005.vcf

    mkdir -p "$SNPDIR"/alt08ref002
    $VARIF -vcf <(gunzip -c "$SNPDIR"/output.vcf.gz) -gff "$GFF" -fasta "$GENOME" \
    --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.8 \
    --ratio-no-alt 0.02 --csv "$SNPDIR"/alt08ref002/filtered-SNPs-sINDELs-08002-best.csv \
    --filteredvcf "$SNPDIR"/alt08ref002/filtered-SNPs-sINDELs-08002-best.vcf

    mkdir -p "$SNPDIR"/alt08ref02
    $VARIF -vcf <(gunzip -c "$SNPDIR"/output.vcf.gz) -gff "$GFF" -fasta "$GENOME" \
    --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.8 \
    --ratio-no-alt 0.2 --csv "$SNPDIR"/alt08ref02/filtered-SNPs-sINDELs-0802-best.csv \
    --filteredvcf "$SNPDIR"/alt08ref02/filtered-SNPs-sINDELs-0802-best.vcf
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
    ./DELLY-caller.sh -v $maxthreads \
        -s "$(ls -1 $BAMFILES | grep -v -E $CONTROLNAME)" -c "$(ls -1 $BAMFILES | grep -E $CONTROLNAME)"
    echo "Done the Other variants step!"
fi
########

##############################--------##############################
