from Bio import SeqIO

input_file = 'out/playground/test.fa'

seqs = []
with open(input_file, 'r') as f:
	seqs = [s for s in SeqIO.parse(f, 'fasta')]

for s in seqs:
	s.seq = s.seq.translate()

with open('out/playground/amino_acids.fa', 'w+') as f:
	SeqIO.write(seqs, f, 'fasta')
