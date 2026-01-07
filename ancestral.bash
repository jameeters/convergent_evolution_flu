OUT_DIR='out/ancestral'

if [[ ! -d "$OUT_DIR" ]]; then
	mkdir "$OUT_DIR"
fi

treetime ancestral \
--tree "out/treetime/timetree.nexus" \
--aln "out/aligned.fa" \
--outdir "$OUT_DIR" \
--report-ambiguous