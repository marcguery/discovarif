#!/bin/bash
#Full pipeline (Wed 30 Sep 11:12:18 CEST 2020)

binversion="0.0.1"

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i doTes=0
declare -i doQua=0
declare -i doMap=0
declare -i doSnp=0
declare -i doCnv=0
declare -i doOth=0
declare -i samplesperrun=1
declare -i threadspersample=1
declare -i dry=0
while getopts ":bqmscok:n:u:h" o; do
    case "${o}" in
        b) # Launch test step.
            doTes=1
            dry=1
            ;;
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
        k) # Use this configuration file
            configfile="$OPTARG"
            ;;
        n) # Process this number of samples in parallel
            samplesperrun=$((OPTARG))
            ;;
        u) # Launch this number of processes per sample
            threadspersample=$((OPTARG))
            ;;
        h | *) # Show help.
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if [[ $SLURM_JOBID =~ ^[0-9]+$ ]] ; then
    LOC=($(echo "$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')"))
    LOC=$(dirname "${LOC[0]}")
    [ $SLURM_CPUS_ON_NODE -lt $((threadspersample*samplesperrun)) ] && {
        echo "Requested $((threadspersample*samplesperrun))"\
        "cores in total but allocated only $SLURM_CPUS_ON_NODE,"\
        "check your slurm configuration"; exit 1; }
else
    LOC=$(dirname "$(realpath $0)")
fi

echo "Discovarif will run with $((threadspersample*samplesperrun)) in total at maximum"

if [ -z "$configfile" -o ! -f "$configfile" ];then
    echo "Selecting the config template file which should only be used for testing"
    echo "Use -k option to provide your own config file"
    configfile="$LOC"/../config/config-template.sh
fi

source "$configfile"

[ ! "$binversion" == "$configversion" ] && { 
    echo "Version of the config file ($configversion)"\
    " does not match the one this pipeline ($binversion)"
    exit 1
 }

[ $samplesperrun -gt $SAMPLENUM ] && { samplesperrun=$SAMPLENUM; }

echo "Samples are: ${SAMPLES[@]}"
##############################-------##############################

##############################PIPELINE##############################

####Test####
if [ $doTes -eq 1 ];then
    echo "Testing each tool..."
    tools=("$FASTQC --version" "$TRIMMOMATIC -version" \
    $BWA "$SAMTOOLS --version" \
    "$BCFTOOLS --version" "$VARIF -h" \
    "$BEDTOOLS --version" "$DELLY --version" "$BEDGRAPHTOBIGWIG" \
    "$PICARD" "$GATK --version")

    for ((i=0;i<${#tools[@]};i++));do
        ${tools[$i]} &> /dev/null
        [ "$?" -eq 127 ] && \
            { echo "Command ${tools[$i]} not found"; echo "Check the executable path in the config file"; }
    done

    echo "Testing raw data files"
    [ -d $DATADIR ] || echo "The folder DATADIR does not exist"
    [ -f $GENOME ] || echo "The file GENOME does not exist"
    [ -f $INFOGENOME ] || echo "The file INFOGENOME does not exist"
    [ -f $GFF ] || echo "The file GFF does not exist"
    [ -f $SAMPLEFILE ] || echo "The file SAMPLEFILE does not exist"
    [ -d $READDIR ] || echo "The folder READDIR does not exist"
    [ -f $CLIPS ] || echo "The file CLIPS does not exist"
    [ -d $OUTDIR ] || echo "The folder OUTDIR does not exist"
    exit 0
fi
########

####Quality with each sample####
if [ $doQua -eq 1 ];then
    echo "Doing the Quality step for each sample..."

    seq -s " " 0 $(($SAMPLENUM-1)) | \
        xargs -d ' ' -n1 -P $samplesperrun bash -c \
        "$LOC"'/../src/mapper-caller.sh -q -s $1 -t '$threadspersample bash
    wait
    echo "Done the Quality step!"
fi
########

####Alignment with each sample####
if [ $doMap -eq 1 ];then
    echo "Doing the Mapping step for each sample..."
    $BWA index "$GENOME"

    seq -s " " 0 $(($SAMPLENUM-1)) | \
        xargs -d ' ' -n1 -P $samplesperrun bash -c \
        "$LOC"'/../src/mapper-caller.sh -m -s $1 -t '$threadspersample bash
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

    seq -s " " 0 $(($SAMPLENUM-1)) | \
        xargs -d ' ' -n1 -P $samplesperrun bash -c \
        "$LOC"'/../src/mapper-caller.sh -v -s $1 -t '$threadspersample bash
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
    
    altlist=("0.8" "0.8" "0.6" "0.51")
    reflist=("0.2" "0.02" "0.4" "0.05")
    for ((i=0;i<$#{altlist[@]};i++));do
        alt=${altlist[$i]}
        ref=${reflist[$i]}
        mkdir -p "$SNPDIR"/alt${alt}ref${ref}

        $VARIF -vcf <(gunzip -c "$SNPDIR"/variants-filtered.vcf.gz) -gff "$GFF" -fasta "$GENOME" \
        --best-variants --all-regions --no-show --depth 5  \
        --ratio-alt ${alt} --ratio-no-alt ${ref} \
        --filtered-csv "$SNPDIR"/alt${alt}ref${ref}/filtered-SNPs-sINDELs.csv \
        --filtered-vcf "$SNPDIR"/alt${alt}ref${ref}/filtered-SNPs-sINDELs.vcf
    done

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
    seq -s " " 0 $(($SAMPLENUM-1)) | \
        xargs -d ' ' -n1 -P $samplesperrun bash -c \
        "$LOC"'/../src/get-coverage.sh -s $1' bash
    wait
    
    ##Getting filtered CNVs
    mkdir -p "$CNVDIR"/summary
    CONTROLSAMPLES=($(awk '$4=="control" { print $1 }' $SAMPLEFILE))
    "$LOC"/../src/CNV-caller.R -indir="$CNVDIR" -incoveragepattern="-perbasecds-core.coverage" \
        -outdir="$CNVDIR"/summary -outcovpattern="-core-cov.tsv" \
        -outsummary="CNV_withoutAPI-MT.csv" -controlsamples="$(echo "${CONTROLSAMPLES[@]}")" \
        -ratiotumor="0.5" -ratiocontrol="0.5"
    
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

    "$LOC"/../src/DELLY-caller.sh -t $((threadspersample*samplesperrun)) -g $GENOME -b "$OUTDIR"/delly-samples.tsv \
        -s "$(echo ${TUMORBAMFILES[@]})" \
        -c "$(echo ${CONTROLBAMFILES[@]})"
    echo "Done the Other variants step!"
fi
########

##############################--------##############################
