#!/bin/bash
#Variant caller pipeline (Mon 28 Sep 11:39:36 CEST 2020)

##############################PARAMETERS##############################
#I had to set those here cause bash would not allow to export array
SAMPLES=($(cut -f1 $SAMPLEFILE | tail -n+2))
READSNAME1=($(cut -f2 $SAMPLEFILE | tail -n+2))
READSNAME2=($(cut -f3 $SAMPLEFILE | tail -n+2))
READS1=($(echo "$(printf $READDIR/'%s\n' "${READSNAME1[@]}")"))
READS2=($(echo "$(printf $READDIR/'%s\n' "${READSNAME2[@]}")"))
#Infer trimmed reads and BAM files
TRIMS1=($(echo "$(printf $TRIMDIR/'%s'_R1$TRIMEXT'\n' "${SAMPLES[@]}")"))
TRIMS2=($(echo "$(printf $TRIMDIR/'%s'_R2$TRIMEXT'\n' "${SAMPLES[@]}")"))
BAMFILES=($(echo "$(printf $BAMBAIDIR/'%s'$BAMEXT'\n' "${SAMPLES[@]}")"))
##############################----------##############################

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i startSample=0
declare -i doQual=0
declare -i doMapp=0
declare -i doVari=0
declare -i threads=1
declare -i memory=1000000000
while getopts ":qmvhs:e:t:g:" o; do
    case "${o}" in
        q) # Launch quality step.
            doQual=1
            ;;
        m) # Launch mapping step.
            doMapp=1
            ;;
        v) # Launch variant calling step (SNP/INDELs)
            doVari=1
            ;;
        s) # Index of the first sample to process.
            startSample=$((OPTARG))
            ;;
        e) # Index of the last sample to process.
            endSample=$((OPTARG))
            ;;
        t) # Use this number of threads for each sample
            threads=$((OPTARG))
            ;;
        g) # Use this memory for each sample
            memory=$((OPTARG))
            ;;
        h | *) # Show help.
            usage
            ;;
    esac
done
shift $((OPTIND-1))
[ -z $endSample ] && { endSample=$(($startSample+1)); }
[ $startSample -ge $endSample -o $endSample -gt $SAMPLENUM ] \
&& { 
    echo "The indices $startSample:$endSample  provided do not match the $SAMPLENUM samples."\
    " Indices should be between 0 and $SAMPLENUM,"\
    " like 0:$SAMPLENUM or 1:$(($SAMPLENUM-1))"
    exit 1; }

bitsToHumanReadable() {
    local i=${1:-0} s=0 S=("" "K" "M" "G" "T")
    while ((i > 1000 && s < ${#S[@]}-1)); do
        i=$((i / 1000))
        s=$((s + 1))
    done
    echo "$i${S[$s]}"
}

memory_format=$(bitsToHumanReadable $memory)
echo "Processing sample indices $startSample to $(($endSample-1))"
##############################-------##############################

##############################PIPELINE##############################
##Data quality with FASTQC/BBDUK
if [ $doQual -eq 1 ];then
    echo "QUALITY step"
    mkdir -p "$QUALDIR"
    mkdir -p "$TRIMDIR"

    for ((i = $startSample ; i < $endSample ; i++ ));do
        r1=${READS1[$i]}
        r2=${READS2[$i]}
        [ ! -f "$r1" -o ! -f "$r2" ] && \
        { echo "At least one read file of sample ${SAMPLES[$i]} does not exist"; continue; }
        r1name=${SAMPLES[$i]}_R1
        r2name=${SAMPLES[$i]}_R2
        echo "Processing sample ${SAMPLES[$i]}..."

        ##Quality raw reads
        $FASTQC -o "$QUALDIR" -t $threads "$r1" "$r2"

        ##Read trimming
        $BBDUK -Xmx"$memory_format" threads=$threads \
            in="$r1" in2="$r2" ref=$CLIPS \
            out="$TRIMDIR/$r1name".paired.gz out2="$TRIMDIR/$r2name".paired.gz \
            outm="$TRIMDIR/$r1name".unpaired.gz outm2="$TRIMDIR/$r2name".unpaired.gz \
            stats="$TRIMDIR/bbduk-stats.txt" \
            k=20 hdist=1 ktrim=r mink=10 qtrim=rl trimq=10 \
            minlength=30 minavgquality=20 minoverlap=14 tbo tpe

        ##Quality trimmed reads
        $FASTQC -o "$QUALDIR" -t $threads \
        "$TRIMDIR/$r1name".paired.gz "$TRIMDIR/$r1name".unpaired.gz \
        "$TRIMDIR/$r2name".paired.gz "$TRIMDIR/$r2name".unpaired.gz
    done
fi

##Mapping with BWA
if [ $doMapp -eq 1 ];then
    echo "MAPPING step"
    mkdir -p "$BAMBAIDIR"
    mkdir -p "$BAMBAIDIR/metrics"
    
    for ((i = $startSample ; i < $endSample ; i++ ));do
        r1=${TRIMS1[$i]}
        r2=${TRIMS2[$i]}
        if [ ! -f "$r1" -o ! -f "$r2" ];then
            echo "Cound not find files $r1 and/or $r2"
            echo "At least one trimmed read file of sample ${SAMPLES[$i]} does not exist, trying with read files provided in $SAMPLEFILE..."
            r1=${READS1[$i]}
            r2=${READS2[$i]}
            if [ ! -f "$r1" -o ! -f "$r2" ];then
                echo "At least one untrimmed read file of sample ${SAMPLES[$i]} does not exist"
                continue
            else
                echo "Found reads from $SAMPLEFILE, beware that using raw reads can lead to bad quality mapping"
            fi
        fi
        samplename=${SAMPLES[$i]}

        echo "Processing sample ${SAMPLES[$i]}..."
        sam="$samplename".sam
        $BWA mem -t $threads "$GENOME" "$r1" "$r2" | \
        $SAMTOOLS addreplacerg -@ $(($threads-1)) -r "ID:${SAMPLES[$i]}" -r "SM:$samplename" - > "$TMPSAM/$sam"

        ##Converting SAM to BAM
        bam="$samplename".bam
        $SAMTOOLS view -@ $(($threads-1)) -S -b "$TMPSAM/$sam" > "$TMPSAM/$bam"

        ##Sorting BAM
        bamsorted="$samplename".sorted.bam
        $SAMTOOLS sort -@ $(($threads-1)) "$TMPSAM/$bam" -o "$TMPSAM/$bamsorted"

        ##Removing duplicates
        bamdedupl="$samplename"$BAMEXT
        java -Xmx"$memory_format" -jar $PICARD_jar MarkDuplicates I="$TMPSAM/$bamsorted" \
            O="$BAMBAIDIR/$bamdedupl" \
            M="$BAMBAIDIR/metrics/$samplename.dd.stats" \
            REMOVE_DUPLICATES=true

        ##Indexing sorted deduplicated BAM
        $SAMTOOLS index -@ $(($threads-1)) "$BAMBAIDIR/$bamdedupl"

        $ALFRED qc -r "$GENOME" -j "$BAMBAIDIR/metrics/${samplename}-qc.json.gz" \
            -o "$BAMBAIDIR/metrics/${samplename}-qc.tsv.gz" "$BAMBAIDIR/$bamdedupl"

        ##Removing tmp files
        [ -s "$BAMBAIDIR/$bamdedupl".bai -a -s "$TMPSAM/$bamsorted" ] && \
            rm "$TMPSAM/$sam" "$TMPSAM/$bam" || \
            echo "Could not remove tmp files because file $BAMBAIDIR/$bamdedupl.bai or $TMPSAM/$bamsorted don't exist or are empty"
    done
fi

if [ $doVari -eq 1 ];then
    #In process of integration
    for ((i = $startSample ; i < $endSample ; i++ ));do
        sample=${BAMFILES[$i]}
        [ ! -f "$sample" ] && \
        { echo "Sample $sample was not found in bam folder $BAMBAIDIR"; continue; }
        samplename=${SAMPLES[$i]}
        
        bamdedupl="$samplename"$BAMEXT
        genomepartsdir="$TMPGVCF/genomeparts"
        gvcf="$samplename".g.vcf.gz

        if [ -f $GVCFDIR/$gvcf -a -f $GVCFDIR/$gvcf.tbi ];then
            echo "$GVCFDIR/$gvcf and $GVCFDIR/$gvcf.tbi already exist, skipping variant calling for sample $samplename"
            continue
        fi

        parts=$((threads/4))

        if [ $parts -lt 1 ];then
            parts=1
        fi
        if [ ! -d "$genomepartsdir/$parts" ];then
            echo "Splitting genome into $parts parts"
            mkdir -p "$genomepartsdir/$parts"
            $SEQKIT split2 -j $threads -p $parts $GENOME -O "$genomepartsdir/$parts"
            for genomepart in $(ls -1 "$genomepartsdir/$parts/"*.fasta);do
                $SEQKIT fx2tab -j $threads -ni "$genomepart" > "$genomepart".list
                rm "$genomepart"
            done
        fi

        partnumber=0
        for genomepart in $(ls -1 "$genomepartsdir/$parts/"*.list);do
            partnumber=$((partnumber+1))

            partgvcf="$samplename"-$partnumber.g.vcf.gz
            partgatkbam="$samplename"-$partnumber.gatk.bam

            { $GATK --java-options "-Xmx${memory_format}" HaplotypeCaller \
                --native-pair-hmm-threads 4 \
                -R "$GENOME" \
                -L "$genomepart" \
                -I "$BAMBAIDIR/$bamdedupl" \
                -O "$TMPGVCF/$partgvcf" \
                -ERC GVCF \
                --sample-ploidy $(($PLOIDY)) \
                -bamout "$BAMBAIDIR/$partgatkbam"; } &
        done
        wait

        java -Xmx"$memory_format" -jar $PICARD_jar MergeVcfs \
            I=<(ls -1 $TMPGVCF/"$samplename"-*.g.vcf.gz) \
            O="$TMPGVCF/$gvcf"
        
        if [ -f $TMPGVCF/$gvcf -a -f "$TMPGVCF/$gvcf.tbi" ];then
            echo "Variant calling complete for sample $samplename, saving in $GVCFDIR"
            cp "$TMPGVCF/$gvcf" "$GVCFDIR/$gvcf"
            cp "$TMPGVCF/$gvcf.tbi" "$GVCFDIR/$gvcf.tbi"
        else
            echo "Failed to call variants for sample $samplename; missing $TMPGVCF/$gvcf or $TMPGVCF/$gvcf.tbi"
        fi
    done

fi

##############################--------##############################
