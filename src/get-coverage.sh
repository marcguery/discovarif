#!/bin/bash
#Coverage pipeline (Wed 30 Sep 11:12:18 CEST 2020)

##############################PARAMETERS##############################
#I had to set those here cause bash would not allow to export array
SAMPLES=($(cut -f1 $SAMPLEFILE | tail -n+2))
BAMFILES=($(echo "$(printf $BAMBAIDIR/'%s'$BAMEXT'\n' "${SAMPLES[@]}")"))
##############################----------##############################

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i startSample=0
while getopts ":hs:e:" o; do
    case "${o}" in
        s) # Index of the first sample to process.
            startSample=$((OPTARG))
            ;;
        e) # Index of the last sample to process.
            endSample=$((OPTARG))
            ;;
        h | *) # Show help.
            usage
            ;;
    esac
done
shift $((OPTIND-1))
[ -z $endSample ] && { endSample=$(($startSample+1)); }
[ $startSample -ge $endSample ] || [ $endSample -gt $SAMPLENUM ] \
&& { 
    echo "The indices $startSample:$endSample  provided do not match the $SAMPLENUM samples."\
    " Indices should be between 0 and $SAMPLENUM,"\
    " like 0:$SAMPLENUM or 1:$(($SAMPLENUM-1))"
    exit 1; }

echo "Processing sample indices $startSample to $(($endSample-1))"
##############################-------##############################

##############################PIPELINE##############################
for ((i = $startSample ; i < $endSample ; i++ ));do
    bamsorted=${BAMFILES[$i]}
    [ ! -f "$bamsorted" ] && \
    { echo "Sample $bamsorted was not found in bam folder $BAMBAIDIR"; continue; }
    samplename=$(cut -d"." -f1 <(basename "$bamsorted"))

    echo "Processing sample ${SAMPLES[$i]}..."
    if [ ! -f "$CNVDIR/$samplename"-cds-core.coverage ];then
        echo "Creating CDS coverage file..."
        $BEDTOOLS coverage -a "$CNVDIR"/view/cds.bed -b "$bamsorted" > "$CNVDIR/$samplename"-cds.coverage
        $BEDTOOLS coverage -a "$CNVDIR"/view/cds-core.bed -b "$bamsorted" > "$CNVDIR/$samplename"-cds-core.coverage
    fi
    
    if [ ! -f "$CNVDIR"/"$samplename"-perbasecds-core.coverage ];then
        echo "Creating perbase CDS coverage file..."
        $BEDTOOLS coverage -a "$CNVDIR"/view/cds.bed -b "$bamsorted" -d > "$CNVDIR/$samplename-"perbasecds.coverage
        $BEDTOOLS coverage -a "$CNVDIR"/view/cds-core.bed -b "$bamsorted" -d > "$CNVDIR/$samplename"-perbasecds-core.coverage
    fi
    
    if [ ! -f "$CNVDIR"/view/"$samplename"-perbase.bed ] && [ ! -f "$CNVDIR"/view/"$samplename"-perbase.bg ];then
        echo "Creating whole perbase file..."
        $BEDTOOLS genomecov -ibam "$bamsorted" -g "$CNVDIR"/chrom.sizes -d > "$CNVDIR"/view/"$samplename"-perbase.bed
    fi

    if [ ! -f "$CNVDIR"/view/"$samplename"-perbase.bg ];then
        echo "Creating whole perbase BG file..."
        awk '{start=$2-1} { printf("%s\t%s\t%i\t%s\n",$1,start,$2,$3); }' "$CNVDIR"/view/"$samplename"-perbase.bed > "$CNVDIR"/view/"$samplename"-perbase.bg
    fi

    if [ ! -f "$CNVDIR"/view/"$samplename"-perbase.bw ];then
        echo "Creating whole perbase BW file..."
        $BEDGRAPHTOBIGWIG "$CNVDIR"/view/"$samplename"-perbase.bg "$CNVDIR"/chrom.sizes "$CNVDIR"/view/"$samplename"-perbase.bw
    fi

    if [ ! -f "$CNVDIR"/view/"$samplename"-perbasecds.bed ];then
        echo "Creating CDS perbase BED file..."
        awk 'FNR==NR{a[$1" "$2]++}FNR!=NR && a[$1" "$2]{print}' "$CNVDIR"/view/cds-perbase.bed "$CNVDIR"/view/"$samplename"-perbase.bed > "$CNVDIR"/view/"$samplename"-perbasecds.bed
    fi
    
    if [ ! -f "$CNVDIR"/view/"$samplename"-perbaseinter.bed ];then
        echo "Creating intergenic perbase BED file..."
        awk 'FNR==NR{a[$1" "$2]++}FNR!=NR && !a[$1" "$2]{print}' "$CNVDIR"/view/cds-perbase.bed "$CNVDIR"/view/"$samplename"-perbase.bed > "$CNVDIR"/view/"$samplename"-perbaseinter.bed
    fi

    if [ ! -f "$CNVDIR"/view/"$samplename"-perbasecds.bg ];then
        echo "Creating CDS perbase BG file..."
        awk '{start=$2-1} { printf("%s\t%s\t%i\t%s\n",$1,start,$2,$3); }' "$CNVDIR"/view/"$samplename"-perbasecds.bed > "$CNVDIR"/view/"$samplename"-perbasecds.bg
    fi

    if [ ! -f "$CNVDIR"/view/"$samplename"-perbasecds.bw ];then
        echo "Creating CDS perbase BW file..."
        $BEDGRAPHTOBIGWIG "$CNVDIR"/view/"$samplename"-perbasecds.bg "$CNVDIR"/chrom.sizes "$CNVDIR"/view/"$samplename"-perbasecds.bw
    fi
    
    if [ ! -f "$CNVDIR"/view/"$samplename"-perbaseinter.bg ];then
        echo "Creating intergenic perbase BG file..."
        awk '{start=$2-1} { printf("%s\t%s\t%i\t%s\n",$1,start,$2,$3); }' "$CNVDIR"/view/"$samplename"-perbaseinter.bed > "$CNVDIR"/view/"$samplename"-perbaseinter.bg
    fi
    
    if [ -f "$CNVDIR"/view/"$samplename"-perbase.bed ] && \
        [ -f "$CNVDIR"/view/"$samplename"-perbasecds.bed ] && \
        [ -f "$CNVDIR"/view/"$samplename"-perbaseinter.bed ];then
        echo "Removing unecessary files..."
        rm "$CNVDIR"/view/"$samplename"-perbase.bed
    fi
    echo "Done with sample ${SAMPLES[$i]}!"
done
##############################--------##############################
