#! /bin/bash

OUT_DIR='out/treetime'

if [[ ! -d "$OUT_DIR" ]]; then
	mkdir "$OUT_DIR"
fi

treetime \
--tree "out/genomic_tree.treefile" \
--aln "out/aligned.fa" \
--dates "out/sample_metadata.tsv" \
--name-column "accession" \
--date-column "stardate" \
--outdir "$OUT_DIR" \
--stochastic-resolve


