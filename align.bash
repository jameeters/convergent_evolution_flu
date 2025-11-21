#! /bin/bash



IN_DIR='/Users/james/Documents/avian-influenza/fasta'
N_SAMPLES=1000

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
| head -n "$N_SAMPLES" \
| xargs cat > "$CONCAT_FA"

# align
ALIGNED_FA="$OUT_DIR/aligned.fa"
mafft --auto "$CONCAT_FA" > "$ALIGNED_FA"