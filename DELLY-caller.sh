#!/bin/bash
#Discover DELLY variants (Wed 30 Sep 18:20:21 CEST 2020)

##############################PARAMETERS##############################
CTRL=$(ls -1 $BAMFILES | grep "$CONTROLNAME")
TUMOR=$(ls -1 $BAMFILES | grep -v "$CONTROLNAME")
##############################----------##############################

##############################PIPELINE##############################

echo "Set $CTRL as controls"
echo "Dealing with $TUMOR..."
variants=("DEL" "DUP" "INV" "TRA" "INS")
if [ ! -f "$DELLYDIR"/MGX-delly-final.vcf ];then
	echo "Variant calling..."
	for variant in ${variants[@]};do
		echo "Doing $variant call..."
        $DELLY call -t $variant -g "$GENOME" "$TUMOR" "$CTRL" \
        -o "$DELLYDIR"/MGX-delly-$variant.bcf && \
		$DELLY filter -t $variant -p -f somatic -m 0 -a 0.5 --ratiogeno 0 -v 5 -c 0.5 \
        "$DELLYDIR"/MGX-delly-$variant.bcf -s $DELLYSAMPLES \
        -o "$DELLYDIR"/MGX-delly-$variant-filter.bcf &
	done
	wait

	echo "Merging files..."
	$BCFTOOLS view $(ls -1 "$DELLYDIR"/*-filter.bcf | head -n1) | grep "#" > "$DELLYDIR"/MGX-delly-final.vcf
	for f in $(ls -1 "$DELLYDIR"/*-filter.bcf);do
		$BCFTOOLS view "$f" | grep -v "#" >> "$DELLYDIR"/MGX-delly-final.vcf
	done
fi

##############################--------##############################
