#!/bin/bash
#Coverage pipeline (Wed 30 Sep 11:12:18 CEST 2020)

##############################PARAMETERS##############################
#I had to set those here cause bash would not allow to export array
SAMPLES=($(ls -1 $READS | grep -o -E $SAMPEXP | uniq))
##############################----------##############################

##############################OPTIONS##############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i startSample=0
declare -i endSample=$SAMPLENUM
while getopts ":hs:" o; do
    case "${o}" in
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

echo "Coverage with sample indices $startSample to $endSample"
##############################-------##############################

##############################PIPELINE##############################
for ((i = $startSample ; i < $endSample ; i++ ));do
    bamsorted=$(grep "${SAMPLES[$i]}" <(ls -1 $BAMFILES))
	[ -z "$bamsorted" ] && \
    { echo "Sample ${SAMPLES[$i]} does not match $BAMFILES"; continue; }
    samplename=$(cut -d"." -f1 <(basename "$bamsorted"))

	echo "Processing sample ${SAMPLES[$i]}..."
	if [ ! -f "$CNVDIR/$samplename"-cds-core.coverage ];then
		echo "Creating CDS coverage file..."
		$BEDTOOLS coverage -a "$CNVDIR"/bed/cds.bed -b "$bamsorted" > "$CNVDIR/$samplename"-cds.coverage
		$BEDTOOLS coverage -a "$CNVDIR"/bed/cds-core.bed -b "$bamsorted" > "$CNVDIR/$samplename"-cds-core.coverage
	fi
	
	if [ ! -f "$CNVDIR"/"$samplename"-perbasecds-core.coverage ];then
		echo "Creating perbase CDS coverage file..."
		$BEDTOOLS coverage -a "$CNVDIR"/bed/cds.bed -b "$bamsorted" -d > "$CNVDIR/$samplename-"perbasecds.coverage
		$BEDTOOLS coverage -a "$CNVDIR"/bed/cds-core.bed -b "$bamsorted" -d > "$CNVDIR/$samplename"-perbasecds-core.coverage
	fi
	
	if [ ! -f "$CNVDIR"/bed/"$samplename"-perbase.bed -a ! -f "$CNVDIR"/bed/"$samplename"-perbasecds.bed -a ! -f "$CNVDIR"/bed/"$samplename"-perbaseinter.bed ];then
		echo "Creating whole perbase file..."
		$BEDTOOLS genomecov -ibam "$bamsorted" -g "$CNVDIR"/chrom.sizes -d > "$CNVDIR"/bed/"$samplename"-perbase.bed
	fi

	if [ ! -f "$CNVDIR"/bed/"$samplename"-perbasecds.bed ];then
		echo "Creating CDS perbase file..."
		awk 'FNR==NR{a[$1" "$2]++}FNR!=NR && a[$1" "$2]{print}' "$CNVDIR"/bed/cds-perbase.bed "$CNVDIR"/bed/"$samplename"-perbase.bed > "$CNVDIR"/bed/"$samplename"-perbasecds.bed
	fi
	
	if [ ! -f "$CNVDIR"/bed/"$samplename"-perbaseinter.bed ];then
		echo "Creating intergenic perbase file..."
		awk 'FNR==NR{a[$1" "$2]++}FNR!=NR && !a[$1" "$2]{print}' "$CNVDIR"/bed/cds-perbase.bed "$CNVDIR"/bed/"$samplename"-perbase.bed > "$CNVDIR"/bed/"$samplename"-perbaseinter.bed
	fi
	
	if [ -f "$CNVDIR"/bed/"$samplename"-perbase.bed -a -f "$CNVDIR"/bed/"$samplename"-perbasecds.bed -a -f "$CNVDIR"/bed/"$samplename"-perbaseinter.bed ];then
		echo "Removing unecessary files..."
		rm "$CNVDIR"/bed/"$samplename"-perbase.bed
	fi
	echo "Done with sample ${SAMPLES[$i]}!"
done
##############################--------##############################
