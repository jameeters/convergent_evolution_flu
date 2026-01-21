#! /bin/bash

INPUT_TREEFILE=$1
ALIGNMENT_FILE=$2
OUT_DIR=$3

if [[ ! -d "$OUT_DIR" ]]; then
	mkdir "$OUT_DIR"
fi

treetime ancestral \
--tree "$INPUT_TREEFILE" \
--aln "$ALIGNMENT_FILE" \
--outdir "$OUT_DIR" \
--report-ambiguous