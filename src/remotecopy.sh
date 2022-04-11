#!/bin/bash
# Download files (Fri 8 Apr 15:54:07 CEST 2022)

#This script will copy your files from a distant server to your current session
# if they are missing from your current session and if you provided a remote address

echo "Copying missing files in $DATADIR from $REMOTEADDRESS:$REMOTEDATADIR"
mkdir -p "$DATADIR"
[ -f $SAMPLEFILE ] && { 
    oldsamplefilename="$(basename "${SAMPLEFILE%.*}".old."${SAMPLEFILE##*.}")"
    echo "Moving current SAMPLEFILE to $(dirname $SAMPLEFILE)/$oldsamplefilename"
    mv $SAMPLEFILE "$(dirname $SAMPLEFILE)/$oldsamplefilename"; }
ssh "$REMOTEADDRESS" [ -d "$REMOTEDATADIR" ] || \
    { echo "$REMOTEDATADIR does not exist in $REMOTEADDRESS"; exit 1; }
rsync -a --ignore-existing --progress "$REMOTEADDRESS":"$REMOTEDATADIR"/ "$DATADIR"/

echo "Copying missing files in $OUTDIR from $REMOTEADDRESS:$REMOTEOUTDIR"
mkdir -p "$OUTDIR"
ssh "$REMOTEADDRESS" [ -d "$REMOTEOUTDIR" ] || \
    { echo "$REMOTEOUTDIR does not exist in $REMOTEADDRESS "; exit 1; }
rsync -a --ignore-existing --progress "$REMOTEADDRESS":"$REMOTEOUTDIR"/ "$OUTDIR"/
