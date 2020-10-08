#!/bin/bash
#Variant caller pipeline (Mon 28 Sep 11:39:36 CEST 2020)

##############################PARAMETERS##############################
#I had to set those here cause bash would not allow to export array
SAMPLES=($(ls -1 $READS | grep -o -E $SAMPEXP | uniq))
PAIREDIDS=( $IDR1 $IDR2 )
##############################----------##############################

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i startSample=0
declare -i endSample=$SAMPLENUM
doQual=0
doMapp=0
doVari=0
while getopts ":qmvhs:" o; do
    case "${o}" in
        q) # Launch quality step.
            doQual=1
            ;;
        m) # Launch mapping step.
            doMapp=1
            ;;
        v) # Launch variant step.
            doVari=1
            ;;
        s) # Indexes of the samples to process.
            startSample=$(cut -d":" -f1 <(echo "${OPTARG}"))
            endSample=$(cut -d":" -f2 <(echo "${OPTARG}"))
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
        r1=$(grep "${SAMPLES[$i]}" <(ls -1 $READS) | grep "${PAIREDIDS[0]}")
        r2=$(grep "${SAMPLES[$i]}" <(ls -1 $READS) | grep "${PAIREDIDS[1]}")
        [ -z "$r1" -o -z "$r2" ] && \
        { echo "Sample ${SAMPLES[$i]} does not match $READS"; continue; }
        r1name=$(cut -d"." -f1 <(basename "$r1"))
        r2name=$(cut -d"." -f1 <(basename "$r2"))
        echo "Processing sample ${SAMPLES[$i]}..."

        ##Quality raw reads
        $FASTQC -o "$QUALDIR" -t 2 "$r1" "$r2"

        ##Read trimming
        $TRIMMOMATIC -trimlog "$TRIMDIR/log.txt" "$r1" "$r2" \
        "$TRIMDIR/$r1name".paired.gz "$TRIMDIR/$r1name".unpaired.gz \
        "$TRIMDIR/$r2name".paired.gz "$TRIMDIR/$r2name".unpaired.gz \
        ILLUMINACLIP:$CLIPS:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

        ##Quality trimmed reads
        $FASTQC -o "$QUALDIR" -t 4 \
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
        r1=$(grep "${SAMPLES[$i]}" <(ls -1 $TRIMREADS) | grep "${PAIREDIDS[0]}")
        r2=$(grep "${SAMPLES[$i]}" <(ls -1 $TRIMREADS) | grep "${PAIREDIDS[1]}")
        [ -z "$r1" -o -z "$r2" ] && \
        { echo "Sample ${SAMPLES[$i]} does not match $TRIMREADS"; continue; }
        samplename=$(cut -d"." -f1 <(basename "$r1") | sed -e 's/'${PAIREDIDS[0]}'//g')
        
        echo "Processing sample ${SAMPLES[$i]}..."
        sam="$samplename".sam
        $BWA index "$GENOME"
        $BWA mem "$GENOME" "$r1" "$r2" > "$SAMDIR/$sam"

        ##Converting SAM to BAM
        bam="$samplename".bam
        $SAMTOOLS view -S -b "$SAMDIR/$sam" > "$BAMBAIDIR/$bam"

        ##Sorting BAM
        bamsorted="$samplename".sorted.bam
        $SAMTOOLS sort "$BAMBAIDIR/$bam" -o "$BAMBAIDIR/$bamsorted"

        ##Indexing sorted BAM
        $SAMTOOLS index "$BAMBAIDIR/$bamsorted"
    done
fi
####UNTESTED####
##Variant calling with bcftools on sorted bams
if [ $doVari -eq 1 ];then
    echo "VARIANT step"
    mkdir -p "$VARIANTDIR"

    for ((i = $startSample ; i < $endSample ; i++ ));do
        bamsorted=$(grep "${SAMPLES[$i]}" <(ls -1 $BAMFILES))
        [ -z "$bamsorted" ] && \
        { echo "Sample ${SAMPLES[$i]} does not match $BAMFILES"; continue; }
        samplename=$(cut -d"." -f1 <(basename $bamsorted))

        echo "Processing sample ${SAMPLES[$i]}..."
        echo $samplename.bcf
        ##Mpileup
        $BCFTOOLS mpileup -f "$GENOME" "$bamsorted" -Ou -d 99999 -a DP,AD,SP | \
        $BCFTOOLS call -m -v -Ob -o "$VARIANTDIR/$samplename".bcf

        ##Filtering Quality superior to 10
        $BCFTOOLS view -I '%QUAL>10' "$VARIANTDIR/$samplename".bcf > \
        "$VARIANTDIR/$samplename"-QUAL10.vcf
    done
fi
####END OF UNTESTED####

##############################--------##############################