#!/bin/bash
#Full pipeline (Wed 30 Sep 11:12:18 CEST 2020)

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i doQua=0
declare -i doMap=0
declare -i doSnp=0
declare -i doCnv=0
declare -i doOth=0
declare -i maxthreads=1
declare -i maxthreadspersample=1
declare -i dry=0
while getopts ":qmscot:u:h" o; do
    case "${o}" in
        q) # Launch quality step.
            doQua=1
            dry=1
            ;;
        m) # Launch mapping step.
            doMap=1
            dry=1
            ;;
        s) # Launch SNP/small INDEL step.
            doSnp=1
            dry=1
            ;;
        c) # Launch CNV step.
            doCnv=1
            dry=1
            ;;
        o) # Launch other variants step.
            doOth=1
            dry=1
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
if [[ $SLURM_JOBID =~ ^[0-9]+$ ]] ; then
    LOC=$(dirname "$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')")
else
    LOC=$(dirname "$(realpath $0)")
fi

source "$LOC"/../config/config.sh

[ $maxthreadspersample -gt $maxthreads ] && \
{ echo "Error: Maximum allowed threads of $maxthreads while $maxthreadspersample requested"; exit 1; }
threadsamples=$(bc <<< $maxthreads/$maxthreadspersample)
[ $threadsamples -gt $SAMPLENUM ] && { threadsamples=$SAMPLENUM; }

echo "Samples are: ${SAMPLES[@]}"
##############################-------##############################

##############################PIPELINE##############################

####Quality with each sample####
if [ $doQua -eq 1 ];then
    echo "Doing the Quality step for each sample..."

    start_index=0
    sequential_per_process=0
    iter=0
    while [ $start_index -le $(($SAMPLENUM-$sequential_per_process)) ];do
        sequential_per_process=$(bc <<< $(($SAMPLENUM-$start_index))/$(($threadsamples-$iter)))
        end_index=$(($start_index+$sequential_per_process))
        echo "New batch: ${SAMPLES[$start_index]} to ${SAMPLES[$(($end_index-1))]}"
        "$LOC"/../src/mapper-caller.sh -q -s $start_index:$end_index -t $maxthreadspersample &
        start_index=$end_index
        ((iter++))
    done
    wait
    echo "Done the Quality step!"
fi
########

####Alignment with each sample####
if [ $doMap -eq 1 ];then
    echo "Doing the Mapping step for each sample..."
    $BWA index "$GENOME"

    start_index=0
    sequential_per_process=0
    iter=0
    while [ $start_index -le $(($SAMPLENUM-$sequential_per_process)) ];do
        sequential_per_process=$(bc <<< $(($SAMPLENUM-$start_index))/$(($threadsamples-$iter)))
        end_index=$(($start_index+$sequential_per_process))
        echo "New batch: ${SAMPLES[$start_index]} to ${SAMPLES[$(($end_index-1))]}"
        "$LOC"/../src/mapper-caller.sh -m -s $start_index:$end_index -t $maxthreadspersample &
        start_index=$end_index
        ((iter++))
    done
    wait
    echo "Done the Mapping step!"
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
    export TMPGVCF=$(mktemp -d "$SNPDIR"/_tmp-gvcf.XXXXXX)
    trap "rm -rf $TMPGVCF" 0 2 3 15

    while [ $start_index -le $(($SAMPLENUM-$sequential_per_process)) ];do
        sequential_per_process=$(bc <<< $(($SAMPLENUM-$start_index))/$(($threadsamples-$iter)))
        end_index=$(($start_index+$sequential_per_process))
        echo "New batch: ${SAMPLES[$start_index]} to ${SAMPLES[$(($end_index-1))]}"
        "$LOC"/../src/mapper-caller.sh -v -s $start_index:$end_index -t $maxthreadspersample &
        start_index=$end_index
        ((iter++))
    done
    wait
    echo "Done the variant step for each sample!"

    $GATK CombineGVCFs \
        -R $GENOME \
        --variant $(sed -e 's/ / --variant /g' <(echo $TMPGVCF/*.g.vcf.gz)) \
        -O $SNPDIR/cohort.g.vcf.gz

    $GATK --java-options "-Xmx4g" GenotypeGVCFs \
        -R $GENOME \
        -V $SNPDIR/cohort.g.vcf.gz \
        -O $SNPDIR/variants.vcf.gz
    
    echo "SNP/INDEL calling terminated"
    echo "Fitering SNP/INDEL..."

    echo "Extracting good quality samples (keep=yes) from the file $SAMPLEFILE..."
    goodsamples=($(awk '$5=="yes" { print $1 }' <(tail -n+2 $SAMPLEFILE)))
    allsamples=($(grep -m 1 "#CHROM" <(gunzip -c "$SNPDIR"/variants.vcf.gz)))
    indices=()
    for sample in ${goodsamples[@]}; do
        samplepresent=1
        i=0
        while [ ! "$sample" == "${allsamples[$i]}" ]; do
            ((i++))
            [ $i -gt ${#allsamples[@]} ] && { echo "Sample $sample is not present in any sample from "$SNPDIR"/variants.vcf.gz, \
                                                you should check the BAM @RG field"; samplepresent=0; break; }
        done
        [ $samplepresent -eq 1 ] && indices+=($((i+1)))
    done
    cut -f 1-9,$(echo ${indices[@]} | sed 's/ /,/g') <(gunzip -c "$SNPDIR"/variants.vcf.gz) | gzip -c > "$SNPDIR"/variants-filtered.vcf.gz
    
    mkdir -p "$SNPDIR"/alt08ref02
    $VARIF -vcf <(gunzip -c "$SNPDIR"/variants-filtered.vcf.gz) -gff "$GFF" -fasta "$GENOME" \
    --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.8 \
    --ratio-no-alt 0.2 --csv "$SNPDIR"/alt08ref02/filtered-SNPs-sINDELs-0802.csv \
    --filteredvcf "$SNPDIR"/alt08ref02/filtered-SNPs-sINDELs-0802.vcf

    ##Filtering variants
    #In renamed vcf files, sample names replace file paths
    mkdir -p "$SNPDIR"/alt06ref04
    $VARIF -vcf <(gunzip -c "$SNPDIR"/variants-filtered.vcf.gz) -gff "$GFF" -fasta "$GENOME" \
    --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.6 \
    --ratio-no-alt 0.4 --csv "$SNPDIR"/alt06ref04/filtered-SNPs-sINDELs-0604.csv \
    --filteredvcf "$SNPDIR"/alt06ref04/filtered-SNPs-sINDELs-0604.vcf

    mkdir -p "$SNPDIR"/alt05ref005
    $VARIF -vcf <(gunzip -c "$SNPDIR"/variants-filtered.vcf.gz) -gff "$GFF" -fasta "$GENOME" \
    --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.55555 \
    --ratio-no-alt 0.05 --csv "$SNPDIR"/alt05ref005/filtered-SNPs-sINDELs-05005.csv \
    --filteredvcf "$SNPDIR"/alt05ref005/filtered-SNPs-sINDELs-05005.vcf

    mkdir -p "$SNPDIR"/alt08ref002
    $VARIF -vcf <(gunzip -c "$SNPDIR"/variants-filtered.vcf.gz) -gff "$GFF" -fasta "$GENOME" \
    --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.8 \
    --ratio-no-alt 0.02 --csv "$SNPDIR"/alt08ref002/filtered-SNPs-sINDELs-08002-best.csv \
    --filteredvcf "$SNPDIR"/alt08ref002/filtered-SNPs-sINDELs-08002-best.vcf

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
    while [ $start_index -le $(($SAMPLENUM-$sequential_per_process)) ];do
        sequential_per_process=$(bc <<< $(($SAMPLENUM-$start_index))/$(($threadsamples-$iter)))
        end_index=$(($start_index+$sequential_per_process))
        echo "New batch: ${SAMPLES[$start_index]} to ${SAMPLES[$(($end_index-1))]}"
        "$LOC"/../src/get-coverage.sh -s $start_index:$end_index &
        start_index=$end_index
        ((iter++))
    done
    wait
    
    ##Getting filtered CNVs
    mkdir -p "$CNVDIR"/summary
    CONTROLSAMPLES=($(awk '$4=="control" { print $1 }' $SAMPLEFILE))
    "$LOC"/../src/CNV-caller.R -indir="$CNVDIR" -incoveragepattern="-perbasecds-core.coverage" \
        -outdir="$CNVDIR"/summary -outcovpattern="-core-cov.tsv" \
        -outsummary="CNV_withoutAPI-MT.csv" -controlsamples="$(echo "${CONTROLSAMPLES[@]}")" \
        -ratiotumor="0.2" -ratiocontrol="1"
    
    echo "Done the CNVs step!"
fi
########

####Other variants####
if [ $doOth -eq 1 ];then
    echo "Doing the Other variants step..."
    mkdir -p "$DELLYDIR"
    cut -f1,4 $SAMPLEFILE | tail -n+2 > "$OUTDIR"/delly-samples.tsv
    CONTROLSAMPLES=($(tail -n+2 $SAMPLEFILE | awk '$4=="control" { print $1 }' $SAMPLEFILE))
    CONTROLBAMFILES=($(echo "$(printf $BAMBAIDIR/'%s'$BAMEXT'\n' "${CONTROLSAMPLES[@]}")"))
    TUMORSAMPLES=($(tail -n+2 $SAMPLEFILE | awk '$4!="control" { print $1 }'))
    TUMORBAMFILES=($(echo "$(printf $BAMBAIDIR/'%s'$BAMEXT'\n' "${TUMORSAMPLES[@]}")"))

    "$LOC"/../src/DELLY-caller.sh -t $maxthreads -g $GENOME -b "$OUTDIR"/delly-samples.tsv \
        -s "$(echo ${TUMORBAMFILES[@]})" \
        -c "$(echo ${CONTROLBAMFILES[@]})"
    echo "Done the Other variants step!"
fi
########

##############################--------##############################
