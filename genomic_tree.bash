#! /bin/bash

IN_FILE="out/ref_aligned_samples_only.fasta"
OUT_PREFIX="out/genomic_tree"

IN_FILE=$1
OUT_PREFIX=$2

iqtree --redo -s "$IN_FILE" --prefix "$OUT_PREFIX"
