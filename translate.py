from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffpandas.gffpandas as gffpd
from pathlib import Path
import subprocess
from Bio.Data.CodonTable import TranslationError
from utils.translation import translate_codonwise
import argparse

parser = argparse.ArgumentParser()

# inputs
parser.add_argument('gff_file')
parser.add_argument('ref_fasta')
parser.add_argument('ancestral_seqs_fasta')

# outputs
parser.add_argument('out_dir')

args = parser.parse_args()

gff_file = args.gff_file
ref_fasta = args.ref_fasta
ancestral_seqs_fasta = args.ancestral_seqs_fasta

out_dir = args.out_dir

Path(out_dir).mkdir(exist_ok=True)



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

with open(f'{out_dir}/coding_seqs.fasta', 'w+') as f:
    SeqIO.write(coding_seqs, f, 'fasta')

aa_seqs = []
for s in coding_seqs:
    translated = translate_codonwise(s)
    aa_seqs.append(translated)

with open(f'{out_dir}/amino_acids.fasta', 'w+') as f:
    SeqIO.write(aa_seqs, f, 'fasta')
