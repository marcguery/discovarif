#!/bin/bash
#Variant caller pipeline (Mon 28 Sep 11:39:36 CEST 2020)

##############################PARAMETERS##############################
#I had to set those here cause bash would not allow to export array
PAIREDIDS=( $IDR1 $IDR2 )
SAMPLES=( $(ls -1 $READS | grep "${PAIREDIDS[0]}" | sort | grep -o -E $SAMPEXP ) )
##############################----------##############################

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i startSample=0
declare -i endSample=$SAMPLENUM
declare -i doQual=0
declare -i doMapp=0
declare -i doVari=0
declare -i threads=1
while getopts ":qmvht:s:" o; do
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
        s) # Indexes of the samples to process.
            startSample=$(cut -d":" -f1 <(echo "${OPTARG}"))
            endSample=$(cut -d":" -f2 <(echo "${OPTARG}"))
            ;;
        t) # Launch this number of processes in parallel
            threads=$((OPTARG))
            ;;
        h | *) # Show help.
            usage
            ;;
    esac
done
shift $((OPTIND-1))
[ $startSample -ge $endSample -o $endSample -gt $SAMPLENUM ] \
&& { 
    echo "The indices provided do not match the $SAMPLENUM samples detected in $READS";
    echo "Indices should be between 0 and $SAMPLENUM, \
like 0:$SAMPLENUM or 1:$(($SAMPLENUM-1))";
    exit 1; }

echo "Variant calling with sample indices $startSample to $endSample"
##############################-------##############################

##############################PIPELINE##############################
##Data quality with FASTQC/TRIMMOMATIC
if [ $doQual -eq 1 ];then
    echo "QUALITY step"
    mkdir -p "$QUALDIR"
    mkdir -p "$TRIMDIR"

    for ((i = $startSample ; i < $endSample ; i++ ));do
        r1=$(ls -1 $READS | grep "${PAIREDIDS[0]}" | sort | head -n $(($i+1)) | tail -n 1 | grep "${SAMPLES[$i]}")
        r2=$(ls -1 $READS | grep "${PAIREDIDS[1]}" | sort | head -n $(($i+1)) | tail -n 1 | grep "${SAMPLES[$i]}")
        [ -z "$r1" -o -z "$r2" ] && \
        { echo "Sample ${SAMPLES[$i]} does not match $READS"; continue; }
        r1name=$(cut -d"." -f1 <(basename "$r1"))
        r2name=$(cut -d"." -f1 <(basename "$r2"))
        echo "Processing sample ${SAMPLES[$i]}..."

        ##Quality raw reads
        $FASTQC -o "$QUALDIR" -t $threads "$r1" "$r2"

        ##Read trimming
        $TRIMMOMATIC -threads $threads -trimlog "$TRIMDIR/log.txt" "$r1" "$r2" \
        "$TRIMDIR/$r1name".paired.gz "$TRIMDIR/$r1name".unpaired.gz \
        "$TRIMDIR/$r2name".paired.gz "$TRIMDIR/$r2name".unpaired.gz \
        ILLUMINACLIP:$CLIPS:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

        ##Quality trimmed reads
        $FASTQC -o "$QUALDIR" -t $threads \
        "$TRIMDIR/$r1name".paired.gz "$TRIMDIR/$r1name".unpaired.gz \
        "$TRIMDIR/$r2name".paired.gz "$TRIMDIR/$r2name".unpaired.gz
    done
fi

##Mapping with BWA
if [ $doMapp -eq 1 ];then
    echo "MAPPING step"
    mkdir -p "$SAMDIR"
    mkdir -p "$BAMBAIDIR"
    
    for ((i = $startSample ; i < $endSample ; i++ ));do
        ##Mapping trimmed paired reads 1 & 2
        r1=$(ls -1 $READS | grep "${PAIREDIDS[0]}" | sort | head -n $(($i+1)) | tail -n 1 | grep "${SAMPLES[$i]}")
        r2=$(ls -1 $READS | grep "${PAIREDIDS[1]}" | sort | head -n $(($i+1)) | tail -n 1 | grep "${SAMPLES[$i]}")
        [ -z "$r1" -o -z "$r2" ] && \
        { echo "Sample ${SAMPLES[$i]} does not match $READS"; continue; }
        samplename=$(cut -d"." -f1 <(basename "$r1") | sed -e 's/'${PAIREDIDS[0]}'//g')

        echo "Processing sample ${SAMPLES[$i]}..."
        sam="$samplename".sam
        $BWA mem -t $threads "$GENOME" "$r1" "$r2" | \
        $SAMTOOLS addreplacerg -@ $(($threads-1)) -r "ID:${SAMPLES[$i]}" -r "SM:$samplename" - > "$SAMDIR/$sam"

        ##Converting SAM to BAM
        bam="$samplename".bam
        $SAMTOOLS view -@ $(($threads-1)) -S -b "$SAMDIR/$sam" > "$BAMBAIDIR/$bam"

        ##Sorting BAM
        bamsorted="$samplename".sorted.bam
        $SAMTOOLS sort -@ $(($threads-1)) "$BAMBAIDIR/$bam" -o "$BAMBAIDIR/$bamsorted"

        ##Removing duplicates
        bamdedupl="$samplename".dd.sorted.bam
        mkdir -p "$BAMBAIDIR/metrics"
        $PICARD MarkDuplicates I="$BAMBAIDIR/$bamsorted" \
            O="$BAMBAIDIR/$bamdedupl" \
            M="$BAMBAIDIR/metrics/$samplename.dd.stats" \
            REMOVE_DUPLICATES=true

        ##Indexing sorted deduplicated BAM
        $SAMTOOLS index -@ $(($threads-1)) "$BAMBAIDIR/$bamdedupl"

        [ $? -eq 0 ] && rm "$BAMBAIDIR/$bamsorted"
    done
fi

if [ $doVari -eq 1 ];then
    #In process of integration
    for ((i = $startSample ; i < $endSample ; i++ ));do
        ##Mapping trimmed paired reads 1 & 2
        sample=$(ls -1 $BAMFILES | sort | head -n $(($i+1)) | tail -n 1 | grep "${SAMPLES[$i]}")
        [ -z "$sample" ] && \
        { echo "Sample ${SAMPLES[$i]} does not match $BAMFILES"; continue; }
        samplename=$(cut -d"." -f1 <(basename "$sample") | sed -e 's/'${PAIREDIDS[0]}'//g')
        
        bamdedupl="$samplename".dd.sorted.bam
        vcf="$samplename".vcf

        echo "$BAMBAIDIR/$bamdedupl"

        continue

        $GATK HaplotypeCaller \
            -R $GENOME \
            -I "$BAMBAIDIR/$bamdedupl" \
            -emit-ref-confidence GVCF \
            --pcr-indel-model NONE \
            --sample-ploidy 1 \
            --max-alternate-alleles 2 \
            --output "$SNPDIR/$vcf"
    done

fi

##############################--------##############################
