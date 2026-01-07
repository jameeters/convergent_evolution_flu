#! /bin/bash

IN_FILE="out/aligned.fa"

iqtree --redo -s "$IN_FILE" --prefix "out/genomic_tree"
