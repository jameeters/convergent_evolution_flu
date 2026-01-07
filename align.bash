#! /bin/bash

SRAS_LIST_FILE="data/sras.txt"

IN_DIR='/Users/james/Documents/avian-influenza/fasta'

OUT_DIR='out'

if [[ ! -d "$OUT_DIR" ]]; then
	mkdir "$OUT_DIR"
fi

# concatenate input sequences
CONCAT_FA="concat.fa"
CONCAT_FA="$OUT_DIR/$CONCAT_FA"

find "$IN_DIR" -name '*HA_cns.fa' \
| sort \
| uniq \
| grep -f "$SRAS_LIST_FILE" \
| xargs cat \
| sed 's/Consensus_//g' \
| sed 's/_HA_cns_threshold_0.75_quality_20//g' \
> "$CONCAT_FA"

# >Consensus_SRR24839058_HA_cns_threshold_0.75_quality_20
# Filter to just the SRA

# align
ALIGNED_FA="$OUT_DIR/aligned.fa"
mafft --auto "$CONCAT_FA" > "$ALIGNED_FA"
