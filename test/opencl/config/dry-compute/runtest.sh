#!/bin/bash

testprog=$1
confdir=$2
outdir=$3

mkdir -p "$outdir"

for conf in "$confdir"/*.yaml; do
    name=${conf%.yaml}
    name=${name##*/}
    echo -n "Running: $name "
    "$testprog" "$conf" "${outdir}/${name}".yaml
done
