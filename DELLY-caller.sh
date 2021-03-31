#!/bin/bash
#Discover DELLY variants (Wed 30 Sep 18:20:21 CEST 2020)

###############################OPTIONS###############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i maxthreads=1
declare -i minsize=0
altaf=0.5
ratiogeno=0
declare -i coverage=5
controlcontamination=2
while getopts ":s:c:t:m:a:r:v:d" o; do
    case "${o}" in
        s) # Files to analyse
            samples="$OPTARG"
            ;;
        c) # Control samples
            controls="$OPTARG"
            ;;
        t) # Number of threads delly will use
            maxthreads=$((OPTARG))
            ;;
        m) # Number of threads delly will use
            minsize=$((OPTARG))
            ;;
        a) # Number of threads delly will use
            altaf=$OPTARG
            ;;
        r) # Number of threads delly will use
            ratiogeno=$OPTARG
            ;;
        v) # Number of threads delly will use
            coverage=$((OPTARG))
            ;;
        d) # Number of threads delly will use
            controlcontamination=$OPTARG
            ;;
    esac
done
shift $((OPTIND-1))
#####################################################################

##############################PARAMETERS##############################
VARIANTS=("DEL" "DUP" "INV" "TRA" "INS")
export OMP_NUM_THREADS=$maxthreads
CTRL=( $controls )
TUMOR=( $samples )
##############################----------##############################

##############################PIPELINE##############################

echo "Set ${CTRL[@]} as controls"
echo "Dealing with ${TUMOR[@]}..."

if [ ! -f "$DELLYDIR"/delly-final.vcf ];then
    echo "Variant calling..."
    for variant in ${VARIANTS[@]};do
        echo "Doing $variant call..."
        $DELLY call -t $variant -g "$GENOME" "${TUMOR[@]}" "${CTRL[@]}" \
        -o "$DELLYDIR"/delly-$variant.bcf && \
        $DELLY filter -t $variant -p -f somatic -m $minsize -a $altaf -r $ratiogeno -v $coverage -c $controlcontamination \
        "$DELLYDIR"/delly-$variant.bcf -s $DELLYSAMPLES \
        -o "$DELLYDIR"/delly-$variant-filter.bcf
    done

    echo "Merging files..."
    $BCFTOOLS view $(ls -1 "$DELLYDIR"/*-filter.bcf | head -n1) | grep "#" > "$DELLYDIR"/delly-final.vcf
    for f in $(ls -1 "$DELLYDIR"/*-filter.bcf);do
        $BCFTOOLS view "$f" | grep -v "#" >> "$DELLYDIR"/delly-final.vcf
    done
fi

##############################--------##############################
