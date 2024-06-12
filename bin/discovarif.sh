#!/bin/bash
#Full pipeline (Wed 30 Sep 11:12:18 CEST 2020)

binversion="0.0.7"
binrealversion="0.0.14"

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i doTes=1
declare -i doQua=0
declare -i doMap=0
declare -i doSnp=0
declare -i doCnv=0
declare -i doOth=0
declare -i samplesperrun=1
declare -i threadspersample=1
declare -i tot_memory=4000000000
declare -i dry=0
while getopts ":qmscok:n:u:g:d:h" o; do
    case "${o}" in
        q) # Launch quality step
            doQua=1
            dry=1
            ;;
        m) # Launch mapping step
            doMap=1
            dry=1
            ;;
        s) # Launch SNP/small INDEL step
            doSnp=1
            dry=1
            ;;
        c) # Launch CNV step
            doCnv=1
            dry=1
            ;;
        o) # Launch other variants step
            doOth=1
            dry=1
            ;;
        k) # Use this configuration file
            configfile="$OPTARG"
            ;;
        n) # Process this number of samples in parallel (default 1)
            samplesperrun=$((OPTARG))
            ;;
        u) # Launch this number of processes per sample (default 1)
            threadspersample=$((OPTARG))
            ;;
        g) # Use this much RAM (default 4G)
            tot_memory=$(($(numfmt --from=si "$OPTARG")))
            ;;
        d) # Path to discovarif 'src' scripts (auto-detect if unset)
            LOC="$OPTARG"
            ;;
        h | *) # Show help
            usage
            ;;
    esac
done
shift $((OPTIND-1))

#######Global configuration#######
echo "discovarif version $binrealversion, global version $binversion"

if [ -z "$configfile" ] || [ ! -f "$configfile" ];then
        echo "You should provide a config file with -k option and fill it accordingly"
        echo "Check in 'config' folder for a template"
        exit 1
fi 

##Run the config file
source "$configfile"
[ ! "$binversion" = "$configversion" ] && { 
    echo "Version of the config file ($configversion)"\
    " is not compatible with the global version ($binversion)"
    exit 1
}
##

##############

#######Check and format cores and RAM#######
#Upper limit of RAM requested
[ $tot_memory -gt 1000000000000000 ] && { tot_memory=1000000000000000; }
#Convert memory to human readable format
bitsToHumanReadable() {
    local i=${1:-0} s=0 S=("" "K" "M" "G" "T")
    while ((i > 1000 && s < ${#S[@]}-1)); do
        i=$((i / 1000))
        s=$((s + 1))
    done
    echo "$i${S[$s]}"
}

tot_memory_format=$(bitsToHumanReadable $tot_memory)

#Halts discovarif if Slurm configuration is invalid
if [[ $SLURM_JOBID =~ ^[0-9]+$ ]];then
    [ $SLURM_CPUS_ON_NODE -lt $((threadspersample*samplesperrun)) ] && {
        echo "Requested $((threadspersample*samplesperrun))"\
        "cores in total but allocated only $SLURM_CPUS_ON_NODE,"\
        "check your slurm configuration"; exit 1; }
fi
##############

#######Find the paths to discovarif scripts#######
declare -i user_input=1
##Set path using calling script when no 'src' directory is set
if [ -z "$LOC" ];then
    user_input=0
    LOC="$(dirname "$0")"
    #Special case when discovarif is called from Slurm
    if [[ $SLURM_JOBID =~ ^[0-9]+$ && ! -f "$LOC"/discovarif.sh ]];then
        LOC=($(echo "$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')"))
        LOC=$(dirname "${LOC[0]}")
    fi
    LOC="$LOC/../src/"
fi
LOC="$(realpath "$LOC")"
##

##Check if src files exist
if [ ! -f "$LOC"/remotecopy.sh ] || [ ! -f "$LOC"/mapper-caller.sh ] || \
   [ ! -f "$LOC"/get-coverage.sh ] || [ ! -f "$LOC"/CNV-caller.R ] || \
   [ ! -f "$LOC"/DELLY-caller.sh ];then
    if [ $user_input -eq 0 ];then
        echo "Auto-detected path '$LOC' does not contain all discovarif scripts"
        echo "Try setting a path using '-d'"
    else
        echo "User-provided path '$LOC' does not contain all discovarif scripts"
        echo "Try unsetting '-d' option to let discovarif auto-detect the path"
    fi
    exit 1
else
    if [ $user_input -eq 0 ];then
        echo "Auto-detected path: $LOC"
    else
        echo "User-provided path: $LOC"
    fi
fi
##############

##############Set up##############

if [ $dry -eq 1 ];then

    ##Copy files if required by config
    if [ ! -z "$REMOTEADDRESS" ];then
        "$LOC"/remotecopy.sh
    fi
    ##

    ##Check integrity
    [ ! -f "$SAMPLEFILE" ] && { echo "Sample file $SAMPLEFILE does not exist"; SAMPLENUM=0; exit 1; }
    newsamplefilename="$(basename "${SAMPLEFILE%.*}".run."${SAMPLEFILE##*.}")"
    head -n1 $SAMPLEFILE > "$(dirname $SAMPLEFILE)/$newsamplefilename"
    tail -n+2 $SAMPLEFILE | awk '$5=="yes" { print $0 }' $SAMPLEFILE >> "$(dirname $SAMPLEFILE)/$newsamplefilename"
    export SAMPLEFILE="$(dirname $SAMPLEFILE)/$newsamplefilename"

    SAMPLES=($(cut -f1 $SAMPLEFILE | tail -n+2))
    #Number of samples
    export SAMPLENUM=$(($(cut -f1 $SAMPLEFILE | tail -n+2 | sort | uniq | wc -l)))

    [ ! $SAMPLENUM -eq ${#SAMPLES[@]} ] &&
    { echo "Found $SAMPLENUM uniquely identified samples but expected ${#SAMPLES[@]}"; \
    echo "The sorted samples with their number of occurences (should all be unique):"; \
    echo "$(cut -f1 $SAMPLEFILE | tail -n+2 | sort | uniq -c)"; exit 1; }

    #Groups
    Groups=($(cut -f4 $SAMPLEFILE | tail -n+2 | sort | uniq))

    ##

    [ $samplesperrun -gt $SAMPLENUM ] && { samplesperrun=$SAMPLENUM; }
    mempersample=$(($tot_memory/$samplesperrun))
    mempersample_format=$(bitsToHumanReadable $mempersample)
    echo "Discovarif will run each sample with with $threadspersample threads and $mempersample_format RAM, totalling to \
    $((threadspersample*samplesperrun)) threads and $tot_memory_format RAM at maximum"

    echo "Samples are: ${SAMPLES[@]}"
fi
############################

##############################-------##############################

##############################PIPELINE##############################

####Test####
if [ $doTes -eq 1 ];then

    echo "Testing each tool..."
    allok=1
    toolok=1
    fileok=1
    tools=("$FASTQC --version" "$BBDUK -version" \
    $BWA "$SAMTOOLS --version" \
    "$BCFTOOLS --version" "$VARIF --version" \
    "$BEDTOOLS --version" "$DELLY --version" "$BEDGRAPHTOBIGWIG" \
    "java -jar $PICARD_jar" "$GATK --version" "$ALFRED --version" \
    "$SEQKIT version")

    for ((i=0;i<${#tools[@]};i++));do
        ${tools[$i]} &> /dev/null
        [ "$?" -eq 127 ] && { toolok=0; allok=0;\
                                echo "Command ${tools[$i]} not found"; \
                                echo "Check the executable path in the config file"; }
    done
    [ $toolok -eq 1 ] && echo "All tools are successfully installed!"

    if [ -z "$REMOTEADDRESS" ] || { [ ! -z "$REMOTEADDRESS" ] && [ $dry -eq 1 ]; };then
        echo "Testing raw data files and folders location..."
        filelist=("$DATADIR" "$GENOME" "$GFF" "$SAMPLEFILE" "$READDIR" "$CLIPS" "$OUTDIR")
        filetype=("d" "f" "f" "f" "d" "f" "d")

        for ((i=0;i<${#filelist[@]};i++));do
            if [ ! -"${filetype[$i]}" "${filelist[$i]}" ];then
                fileok=0
                allok=0
                echo "The file or directory ${filelist[$i]} does not exist"
            fi
        done
    else
        fileok=0
        echo "Files and folders are located on a remote server ($REMOTEADDRESS)"
        echo "Cannot test locations without running one step of the pipeline"
        exit 1
    fi
    [ $fileok -eq 1 ] && echo "All files and folders were found!"
    [ $allok -eq 1 ] && { echo "Config file is properly configured!"; } || \
        { echo "Config file is not properly configured"; exit 1; }
fi
########

####Quality with each sample####
if [ $doQua -eq 1 ];then
    echo "Doing the Quality step for each sample..."

    seq -s " " 0 $(($SAMPLENUM-1)) | \
        xargs -d ' ' -n1 -P $samplesperrun bash -c \
        "$LOC"'/mapper-caller.sh -q -s $1 -t '$threadspersample' -g '$mempersample bash
    wait
    echo "Done the Quality step!"
fi
########

####Alignment with each sample####
if [ $doMap -eq 1 ];then
    echo "Doing the Mapping step for each sample..."
    export TMPSAM=$(mktemp -d "$TMPSAMDIR"/sam.XXXXXX)
    trap 'rm -rf -- "$TMPSAM"' EXIT
    echo "The temporary folder for this mapping run is $TMPSAM"

    $BWA index "$GENOME"

    seq -s " " 0 $(($SAMPLENUM-1)) | \
        xargs -d ' ' -n1 -P $samplesperrun bash -c \
        "$LOC"'/mapper-caller.sh -m -s $1 -t '$threadspersample' -g '$mempersample bash
    wait
    echo "Done the Mapping step!"
fi
########

####SNP/small INDEL filtering####
if [ $doSnp -eq 1 ];then
    echo "Doing the SNPs/INDELs step..."
    mkdir -p "$SNPDIR"
    mkdir -p "$SNPDIR"/_tmp/
    mkdir -p "$SNPDIR"/varif_output/
    mkdir -p "$GVCFDIR"

    $GATK --java-options "-Xmx${tot_memory_format}" CreateSequenceDictionary \
        -R "$GENOME" &>/dev/null
    
    $SAMTOOLS faidx "$GENOME"

    export TMPGVCF=$(mktemp -d "$TMPGVCFDIR"/gvcf.XXXXXX)
    trap 'rm -rf -- "$TMPGVCF"' EXIT
    echo "The temporary folder for this GATK run is $TMPGVCF"

    seq -s " " 0 $(($SAMPLENUM-1)) | \
        xargs -d ' ' -n1 -P $samplesperrun bash -c \
        "$LOC"'/mapper-caller.sh -v -s $1 -t '$threadspersample' -g '$mempersample bash
    wait
    echo "Done the variant step for each sample!"

    $GATK --java-options "-Xmx${tot_memory_format}" CombineGVCFs \
        -R $GENOME \
        --variant $(sed -e 's/ / --variant /g' <(echo $GVCFDIR/*.g.vcf.gz)) \
        -O $SNPDIR/cohort.g.vcf.gz

    $GATK --java-options "-Xmx${tot_memory_format}" GenotypeGVCFs \
        -R $GENOME \
        -V $SNPDIR/cohort.g.vcf.gz \
        -O $SNPDIR/variants.vcf.gz
    
    echo "SNP/INDEL calling terminated"

    #Samples are sorted according to their group identifier (must be numeric)
    echo "Extracting good quality samples (keep=yes) from the file $SAMPLEFILE..."
    goodsamples=($(awk '$5=="yes" { print $1 }' <(tail -n+2 $SAMPLEFILE |  sort -k4 -n)))
    goodsamplesgroups=($(awk '$5=="yes" { print $4 }' <(tail -n+2 $SAMPLEFILE |  sort -k4 -n)))
    allsamples=($(grep -m 1 "#CHROM" <(gunzip -c "$SNPDIR"/variants.vcf.gz) | cut -f10-))
    goodindices=()
    allgoodindices=()
    group=
    controlindices=()
    #Extract samples from VCF file according to their groups
    for ((i=0;i<${#goodsamples[@]};i++));do
        sample=${goodsamples[$i]}
        samplegroup=${goodsamplesgroups[$i]}
        samplepresent=1
        indexvcf=0
        #Iterate over all VCF headers to check if the sample is present
        while [ ! "$sample" = "${allsamples[$indexvcf]}" ]; do
            ((indexvcf++))
            [ $indexvcf -gt ${#allsamples[@]} ] && { 
                echo "Sample $sample is not present in any sample from "$SNPDIR"/variants.vcf.gz,"\
                " you should check the BAM @RG field"
                samplepresent=0
                break; }
        done
        [ $samplepresent -eq 1 ] && {
            #If group is set and is different from the samplegroup (the current one), ...
            if [ ! -z "$group" ] && [ "$group" != "$samplegroup" ];then
                #if there is a control group (must be '0') it is the first one
                # to appear since samples are sorted by their group id
                [ "$group" = "0" ] && {
                    echo "Found control group 0 at columns ${goodindices[@]} of "$SNPDIR"/variants.vcf.gz"
                    controlindices=("${goodindices[@]}")
                }
                #... group samples are extracted with the control samples (if '0' id is present) from the VCF file
                goodindices+=(${controlindices[@]})
                echo "Group $group will be constituted by columns $(echo "${goodindices[@]}" | sed 's/ /,/g') of "$SNPDIR"/variants.vcf.gz"
                cut -f 1-9,$(echo "${goodindices[@]}" | sed 's/ /,/g') <(gunzip -c "$SNPDIR"/variants.vcf.gz) | gzip -c > "$SNPDIR"/variants-filtered-$group.vcf.gz
                goodindices=()
            fi
            goodindices+=($((indexvcf+10)))
            allgoodindices+=($((indexvcf+10)))
            #group is now samplegroup (current group)
            group=$samplegroup
        }
    done
    goodindices+=(${controlindices[@]})
    echo "Group $group will be constituted by columns $(echo "${goodindices[@]}" | sed 's/ /,/g') of "$SNPDIR"/variants.vcf.gz"
    cut -f 1-9,$(echo "${goodindices[@]}" | sed 's/ /,/g') <(gunzip -c "$SNPDIR"/variants.vcf.gz) | gzip -c > "$SNPDIR"/variants-filtered-$group.vcf.gz

    allgoodindices+=(${controlindices[@]})
    echo "All available samples are constituted by columns $(echo "${allgoodindices[@]}" | sed 's/ /,/g') of "$SNPDIR"/variants.vcf.gz"
    cut -f 1-9,$(echo "${allgoodindices[@]}" | sed 's/ /,/g') <(gunzip -c "$SNPDIR"/variants.vcf.gz) | gzip -c > "$SNPDIR"/variants-filtered.vcf.gz

    echo "Filtering SNPs and INDELs with varif"
    gzip -c "$SNPDIR"/variants-filtered.vcf.gz > "$SNPDIR"/variants-filtered.vcf

    varifsamples=$(mktemp varif-samples.XXXXXX.tsv)
    echo "The temporary file for this varif run is $OUTDIR/${varifsamples}"
    awk -v OFS=$'\t' '{print ($4,$1,"0","0","other","0")}' $SAMPLEFILE | tail -n+2 > "$OUTDIR/${varifsamples}"

    $VARIF -vcf "$SNPDIR"/variants-filtered.vcf -gff "$GFF" -fasta "$GENOME" \
            -outfilename "$SNPDIR"/varif_output/filtered-SNPs-sINDELs-alt0.9-ref0.05 \
            --comparison all \
            --best-variants --all-regions --depth 5 \
            --ratio-alt 0.9 --ratio-ref 0.05 --output-vcf \
            --ncores $((threadspersample*samplesperrun)) --chunk-size 5000
    
    $VARIF -vcf "$SNPDIR"/variants-filtered.vcf -gff "$GFF" -fasta "$GENOME" \
            -outfilename "$SNPDIR"/varif_output/filtered-SNPs-sINDELs-alt0.9-ref0.05 \
            --ped "$OUTDIR/${varifsamples}" --comparison families \
            --best-variants --all-regions --depth 5 \
            --ratio-alt 0.9 --ratio-ref 0.05 --output-vcf \
            --ncores $((threadspersample*samplesperrun)) --chunk-size 5000
    
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
        "$LOC"'/get-coverage.sh -s $1' bash
    wait
    
    ##Getting filtered CNVs
    mkdir -p "$CNVDIR"/summary
    CONTROLSAMPLES=($(awk '$4=="control" { print $1 }' $SAMPLEFILE))
    "$LOC"/CNV-caller.R -indir="$CNVDIR" -incoveragepattern="-perbasecds-core.coverage" \
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
    dellysamples=$(mktemp "$OUTDIR/"delly-samples.XXXXXX.tsv)
    echo "The temporary file for this DELLY run is $dellysamples"

    awk '{if ($4==0) ($4 = "control"); else ($4 = "tumor"); print ($1,$4)}' $SAMPLEFILE | tail -n+2 > "${dellysamples}"
    CONTROLSAMPLES=($(tail -n+2 $SAMPLEFILE | awk '$4=="0" { print $1 }' $SAMPLEFILE))
    CONTROLBAMFILES=($(echo "$(printf $BAMBAIDIR/'%s'$BAMEXT'\n' "${CONTROLSAMPLES[@]}")"))
    TUMORSAMPLES=($(tail -n+2 $SAMPLEFILE | awk '$4!="0" { print $1 }'))
    TUMORBAMFILES=($(echo "$(printf $BAMBAIDIR/'%s'$BAMEXT'\n' "${TUMORSAMPLES[@]}")"))
    VARIANTS=("DEL" "DUP" "INS" "INV" "BND")
    if [ $samplesperrun -gt ${#VARIANTS[@]} ];then
        dellythreads=$(((threadspersample*samplesperrun)/${#VARIANTS[@]}))
    else
        dellythreads=$threadspersample
    fi

    echo -n "${VARIANTS[@]}" | \
        xargs -d ' ' -n1 -P $samplesperrun bash -c \
            '$1 -t "$2" -g "$3" -b "$4" -s "$5" -c "$6" -u "$7" ' \
        bash "$LOC"/DELLY-caller.sh $dellythreads "$GENOME" "${dellysamples}" "$(echo ${TUMORBAMFILES[@]})" "$(echo ${CONTROLBAMFILES[@]})"
    wait
    
    echo "Merging files..."
    compgen -G "$DELLYDIR"/*-filter.bcf > /dev/null
    returnerr=$?
    if [ ! $returnerr -eq 0 ];then
        echo "Cannot merge filtered files from DELLY"
        exit $returnerr
    fi
    bcffiles=("$DELLYDIR"/*-filter.bcf)
    $BCFTOOLS view "${bcffiles[0]}" | grep "#" > "$DELLYDIR"/delly-final.vcf
    for f in "$DELLYDIR"/*-filter.bcf;do
        $BCFTOOLS view "$f" | grep -v "#" >> "$DELLYDIR"/delly-final.vcf
    done

    echo "Done the Other variants step!"
fi
########

##############################--------##############################
