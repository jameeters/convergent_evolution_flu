#! /bin/bash


TREETIME_OUT_DIR="out/ancestral"

treetime homoplasy \
--tree "$TREETIME_OUT_DIR/annotated_tree.nexus" \
--aln "out/aligned.fa" 
