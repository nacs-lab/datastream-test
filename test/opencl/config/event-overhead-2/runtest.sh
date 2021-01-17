#!/bin/bash

testprog=$1
outdir=$2
confdir=${3:-$(dirname "${BASH_SOURCE}")/generic}

mkdir -p "$outdir"

for conf in "$confdir"/*.yaml; do
    name=${conf%.yaml}
    name=${name##*/}
    echo -n "Running: $name "
    "$testprog" "$conf" | tee "${outdir}/${name}".yaml
done
