# Pipeline Description

Last updated: 13 Jan 2026

## Steps

### Setup

`data/sras.txt`: must list the SRAs for samples to be analyzed.

### Get sample metadata

`get_sample_metadata.py`

**Inputs:**

- `sras.txt`

**Outputs:**

- `sample_metadata.tsv`

Fetch sample metadata from Muninn (kenny).
Must be run before treetime, as it provides collection dates for the samples.

### Multiple sequence alignment

`align.bash`

**Inputs:**

- `sras.txt`
- Fasta files for each sample

**Outputs:**

- `concat.fa`
- `aligned.fa`

Concatenate all sample sequences into a single file, then align using mafft.

### Genomic tree

`genomic_tree.bash`

**Inputs:**

- `aligned.fa`

**Outputs:**

- `genomic_tree.*`

Constructs a ML phylogeny given the aligned samples

### Time-scaled tree

`time_tree.bash`

**Inputs:**

- `genomic_tree.treefile`
- `aligned.fa`
- `sample_metadata.tsv`

**Outputs:**

- `out/treetime/*`

Given the genomic tree, use TreeTime to produce a time-scaled tree.

### Ancestral sequence reconstruction

`ancestral.bash`

**Inputs:**

- `aligned.fa`
- `treetime/timetree.nexus`

**Outputs:**

- `out/ancestral/*`

Given the time-scaled tree, apply TreeTime again to reconstruct genomes for inferred ancestors.

### Translation to amino acid sequences

`translate.py`

**Inputs:**

- `data/HA.gff3` gff file describing CDS regions of the gene
- `data/HA_reference.fasta` nucleotide sequence of the reference version of the gene, corresponding to the GFF file.
- `ancestral/ancestral_sequences.fasta`

**Outputs:**

- `out/translation/*`
- `ref_aligned_nt.fasta`
- `coding_seqs.fasta`
- `amino_acids.fasta`

Align all sample and ancestral sequences (`ancestral_sequences.fasta`) to a known reference, corresponding to a GFF file describing the coding regions of the gene.
Write result to `ref_aligned_nt.fasta`.

Read GFF file and extract starts and ends for each CDS in the gene.
Subset the ref-aligned sequences to coding regions only (output: `coding_seqs.fasta`).

Translate the coding sequences into amino acid sequences (`amino_acids.fasta`).

### Analyze tree, count mutation occurrences

`read_tree.py`

**Inputs:**

- `ancestral/annotated_tree.nexus`
- `translation/coding_seqs.fasta`
- `translation/amino_acids.fasta`

**Outputs:**

- `/tmp/drawing.svg` Figure of the tree, currently not used for anything.

Read ancestral tree.
Read coding sequences and amino acids files and assign sequences to each node in the tree.
Compare the coding sequence of each node to that of its parent codon-by-codon and note mutations.
Filter out mutations that lead to an X.
Count the synonymous and non-synonymous mutations at each site, and for sites that have both, compute the dn/ds ratio.
Because we are starting with a node's mutations relative to its parent, we are counting the number of times that each mutation has arisen in our tree.

For now the amino acid sequences are only used to help check the results of the codon comparison.
