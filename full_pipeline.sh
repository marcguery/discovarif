#!/bin/bash

##############################PIPELINE##############################
source config.sh

#Quality and alignment with each sample
./variant-caller.sh -qm -s 0:2 &
./variant-caller.sh -qm -s 2:4 &
./variant-caller.sh -qm -s 4:6 &
./variant-caller.sh -qm -s 6:8 &
./variant-caller.sh -qm -s 8:10 &
./variant-caller.sh -qm -s 10:12 &
./variant-caller.sh -qm -s 12:14 &
wait
exit 0

#Doing the mpileup with all the samples at once
$BCFTOOLS mpileup --threads 7 -f $GENOME "$BAMBAIDIR"/*.sorted.bam -Ou -d 99999 -a DP,AD,SP | \
$BCFTOOLS call -m -v -Ob -o "$VARIANTDIR/"all-samples.bcf

$BCFTOOLS view --threads 7 -i '%QUAL>=10' "$VARIANTDIR/"all-samples.bcf > "$VARIANTDIR/"all-samples-Qual10.vcf

##SNP/INDEL filtering

$VARIF -vcf "$VARIANTDIR/"all-samples-Qual10.vcf -gff $GFF -fasta $GENOME \
--no-fixed --best-variants --all-regions --no-show --depth 6 --ratio-alt 0.8 \
--ratio-no-alt 0.2 --csv "$VARIANTDIR/"all-samples-SNPs.csv

##CNV

##Other variants


##############################--------##############################