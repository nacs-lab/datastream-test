#!/bin/bash -e

diffprog=$1
ifile=$2
ofile=$3

ls -l "$ifile"
zstd -d "$ifile" -o "$ifile".tmp
"$diffprog" "$ifile".tmp "$ifile".tmp
zstd -f -19 "$ifile".tmp -o "$ofile"
rm "$ifile".tmp
