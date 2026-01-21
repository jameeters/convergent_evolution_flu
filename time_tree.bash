#! /bin/bash

GENOMIC_TREEFILE=$1
ALIGNMENT_FILE=$2
METADATA_FILE=$3
OUT_DIR=$4

if [[ ! -d "$OUT_DIR" ]]; then
	mkdir "$OUT_DIR"
fi

treetime \
--tree "$GENOMIC_TREEFILE" \
--aln "$ALIGNMENT_FILE" \
--dates "$METADATA_FILE" \
--name-column "accession" \
--date-column "stardate" \
--outdir "$OUT_DIR" \
--stochastic-resolve


