#!/bin/sh

if [ $# -ne 1 ]; then
    echo "Usage: GenICB2.sh <input_dir>" 1>&2
    exit 1
fi

echo "Input dir:  $1"

for ifile in `\find $1 -maxdepth 1 -type f`; do
    ofile="$ifile.bin"
    # echo "$ifile $ofile"
    ./GenICB2 $ifile $ofile
done
