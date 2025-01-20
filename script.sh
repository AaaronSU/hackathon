#!/bin/bash

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <filename> <#simulations> <#runs>"
    exit 1
fi

filename=$1
n_simulations=$2
n_runs=$3

dir="bench"
mkdir -p $dir
mkdir -p ${dir}/all
mkdir -p ${dir}/resume

CC=g++
OFLAGS_ARRAY=("-O0" "-O1" "-O2" "-O3" "-Ofast")
FLAGS="-march=native -funroll-loops" # Pour ARM, x86_64 = -march=native -funroll-loops

# Compilation avec g++

## Sans flags d'optimisation
$CC -O $filename.cxx -o $filename

## Exécution
output="$dir/all/output_ref.txt"
echo "Simulations = $n_simulations, Runs = $n_runs" > "$output"
for i in {1..10}; do
    echo "Run $i" >> $output
    ./$filename $n_simulations $n_runs >> $output
    echo " " >> $output
done

## Avec ciblage de l'architecture + -funroll-loops et =! FLAGS
for flag in ${OFLAGS_ARRAY[@]}; do
    output="$dir/all/$output_${flag//-/}.txt"
    echo "Simulations = $n_simulations, Runs = $n_runs" > "$output"

    $CC $OFLAGS $FLAGS $filename.cxx -o $filename

    ## Exécution
    for i in {1..10}; do 
        echo "Run $i" >> $output
        ./$filename $n_simulations $n_runs >> $output
        echo " " >> $output
    done
done

# Lecture du fichier output.txt, extraction de %Average relative error:
for bench_file in ${dir}/all/*; do
    name=$(basename $bench_file)
    output="${dir}/resume/${name%.*}.txt"
    echo "$n_simulations $n_runs" > $output 
    grep "%Average Relative Error:" "$bench_file" | sed 's/.*%Average Relative Error: //' >> "$output"
done