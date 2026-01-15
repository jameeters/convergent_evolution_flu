from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffpandas.gffpandas as gffpd
from pathlib import Path
import subprocess



gff_file = 'data/HA.gff3'
ref_fasta = 'data/HA_reference.fasta'

ancestral_seqs_fasta = 'out/ancestral/ancestral_sequences.fasta'

Path('out/translation').mkdir(exist_ok=True)


# Read GFF, get slice indices for CDS
cds_data = gffpd.read_gff3(gff_file).filter_feature_of_type(['CDS'])

slices = []
for row in cds_data.df.iterrows():
    cds = row[1]
    start = cds['start']
    end = cds['end']
    # convert to slice coordinates: zero indexed, half open
    slices.append((start-1, end))

# read aligned nt seqs
seqs = []
with open(ancestral_seqs_fasta, 'r') as f:
    seqs = [s for s in SeqIO.parse(f, 'fasta')]

# apply slices based on CDS and translate
coding_seqs = []
for s in seqs:
    coding = Seq('')
    for start, end in slices:
        coding += s.seq[start:end]
    coding_seqs.append(SeqRecord(coding, id=s.id, name=s.name, description=s.description))

with open('out/translation/coding_seqs.fasta', 'w+') as f:
    SeqIO.write(coding_seqs, f, 'fasta')

aa_seqs = []
for s in coding_seqs:
    translated = SeqRecord(s.seq.translate(), id=s.id, name=s.name, description=s.description)
    aa_seqs.append(translated)

with open('out/translation/amino_acids.fasta', 'w+') as f:
    SeqIO.write(aa_seqs, f, 'fasta')
