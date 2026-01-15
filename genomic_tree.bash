#! /bin/bash

IN_FILE="out/ref_aligned_samples_only.fasta"

iqtree --redo -s "$IN_FILE" --prefix "out/genomic_tree"
