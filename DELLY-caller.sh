#!/bin/bash
#Discover DELLY variants (Wed 30 Sep 18:20:21 CEST 2020)
DELLY="${DELLY:-delly}"
BCFTOOLS="${BCFTOOLS:-bcftools}"
DELLYDIR="${DELLYDIR:-.}"
###############################OPTIONS###############################
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }
declare -i maxthreads=1
declare -i minsize=0
altaf=0.5
ratiogeno=0
declare -i coverage=5
controlcontamination=2
while getopts ":s:c:b:g:t:m:a:r:v:d:" o; do
    case "${o}" in
        s) # Files to analyse (required)
            samples="$OPTARG"
            ;;
        c) # Control samples (required)
            controls="$OPTARG"
            ;;
        b) # Tabulated DELLY config file (required)
            dellysamples="$OPTARG"
            ;;
        g) # Genome in a fasta format (required)
            genome="$OPTARG"
            ;;
        t) # Number of threads delly will use
            maxthreads=$((OPTARG))
            ;;
        m) # Minimal size of a variant for it te be kept
            minsize=$((OPTARG))
            ;;
        a) # Alternative allele frequency for a variant to be considered true
            altaf=$OPTARG
            ;;
        r) # Fraction of the samples having enough coverage for a variant to be annotated 'PASS'
            ratiogeno=$OPTARG
            ;;
        v) # Minimal coverage to annotate a sample variant as 'PASS'
            coverage=$((OPTARG))
            ;;
        d) # Maximal fraction of variant reads allowed on the control samples
            controlcontamination=$OPTARG
            ;;
    esac
done
echo "Arguments were: $samples / $controls / $dellysamples / $genome / $maxthreads / $minsize / $altaf / $ratiogeno / $coverage / $controlcontamination"
shift $((OPTIND-1))
[ -z "$controls" -o -z "$samples" -o -z "$dellysamples" -o -z "$genome" ] && { echo "Missing arguments, aborting..."; usage; }
#####################################################################

##############################PARAMETERS##############################
VARIANTS=("DEL" "DUP" "INS" "INV" "BND")
export OMP_NUM_THREADS=$maxthreads
CTRL=( $controls )
TUMOR=( $samples )
declare -i oktocontinue=0

for file in "${TUMOR[@]}"; do
    [ ! -f $file ] && { echo "Sample file $file does not exist"; oktocontinue=1; }
done
for file in "${CTRL[@]}"; do
    [ ! -f $file ] && { echo "Control file $file does not exist"; oktocontinue=1; }
done
[ ! -f $dellysamples ] && { echo "DELLY config file $dellysamples does not exist"; oktocontinue=1; }
[ ! -f $genome ] && { echo "Genome file $genome does not exist"; oktocontinue=1; }

[ "$oktocontinue" -eq 1 ] && exit 1
##############################----------##############################

##############################PIPELINE##############################

echo "Running DELLY with $OMP_NUM_THREADS threads"
echo "Set ${CTRL[@]} as controls"
echo "Dealing with ${TUMOR[@]}..."

echo "Variant calling..."
for variant in ${VARIANTS[@]};do
    if [ ! -f "$DELLYDIR"/delly-$variant.bcf ];then
        echo "Doing $variant call..."
        $DELLY call -t $variant -g "$genome" "${TUMOR[@]}" "${CTRL[@]}" \
           -o "$DELLYDIR"/delly-$variant.bcf
    fi
    if [ ! -f "$DELLYDIR"/delly-$variant-filter.bcf ];then
        echo "Doing variant filtering..."
        $DELLY filter -p -f somatic -m $minsize -a $altaf -r $ratiogeno -v $coverage -c $controlcontamination \
            "$DELLYDIR"/delly-$variant.bcf -s $dellysamples \
            -o "$DELLYDIR"/delly-$variant-filter.bcf
        [ ! -f "$DELLYDIR"/delly-$variant-filter.bcf ] && touch "$DELLYDIR"/delly-$variant-filter.bcf
    fi
done

    echo "Merging files..."
    $BCFTOOLS view $(ls -1 "$DELLYDIR"/*-filter.bcf | head -n1) | grep "#" > "$DELLYDIR"/delly-final.vcf
    for f in $(ls -1 "$DELLYDIR"/*-filter.bcf);do
        $BCFTOOLS view "$f" | grep -v "#" >> "$DELLYDIR"/delly-final.vcf
    done
fi

##############################--------##############################
