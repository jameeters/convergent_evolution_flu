#! /bin/bash

set -euxo pipefail

OUT_DIR='out_bovine_b3.13_pb2'

SRAS_FILE="data/sras_bovine_b3.13.txt"
FASTAS_DIR="/Users/james/Documents/avian-influenza/fasta"
REF_FASTA="data/PB2_reference.fasta"
REGION_NAME="PB2"

CONCATENATED_FASTA="$OUT_DIR/concat.fasta"
SRAS_QC_PASSED_FILE="$OUT_DIR/sras_qc_passed.txt"
ALIGNED_WITH_REF_FASTA="$OUT_DIR/ref_aligned.fasta"
ALIGNED_SAMPLES_ONLY_FASTA="$OUT_DIR/ref_aligned_samples_only.fasta"


python3 align.py \
"$SRAS_FILE" \
"$FASTAS_DIR" \
"$REF_FASTA" \
"$REGION_NAME" \
"$OUT_DIR" \
"$ALIGNED_SAMPLES_ONLY_FASTA"

SAMPLE_METADATA_FILE="$OUT_DIR/sample_metadata.tsv"


python3 get_sample_metadata.py \
"$SRAS_FILE" \
"$SAMPLE_METADATA_FILE"

GENOMIC_TREE_OUTPUT_PREFIX="$OUT_DIR/genomic_tree"

./genomic_tree.bash \
"$ALIGNED_SAMPLES_ONLY_FASTA" \
"$GENOMIC_TREE_OUTPUT_PREFIX"

TREETIME_OUT_DIR="$OUT_DIR/treetime"

./time_tree.bash \
"$GENOMIC_TREE_OUTPUT_PREFIX.treefile" \
"$ALIGNED_SAMPLES_ONLY_FASTA" \
"$SAMPLE_METADATA_FILE" \
"$TREETIME_OUT_DIR"

ANCESTRAL_OUT_DIR="$OUT_DIR/ancestral"

./ancestral.bash \
"$TREETIME_OUT_DIR/timetree.nexus" \
"$ALIGNED_SAMPLES_ONLY_FASTA" \
"$ANCESTRAL_OUT_DIR"


HA_GFF_FILE="data/PB2.gff3"
TRANSLATION_OUT_DIR="$OUT_DIR/translation"

python3 translate.py \
"$HA_GFF_FILE" \
"$REF_FASTA" \
"$ANCESTRAL_OUT_DIR/ancestral_sequences.fasta" \
"$TRANSLATION_OUT_DIR"

python3 read_tree.py \
"$ANCESTRAL_OUT_DIR" \
"$TRANSLATION_OUT_DIR" \
"$OUT_DIR"
