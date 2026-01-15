from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffpandas.gffpandas as gffpd
from pathlib import Path
import subprocess
import glob
import os

# inputs
sras_file = "data/sras.txt"
fastas_dir = '/Users/james/Documents/avian-influenza/fasta'
ambiguous_bases_frac_limit = 0.01
ref_fasta = 'data/HA_reference.fasta'

# outputs todo: remove test
out_dir = 'out'
concatenated_fasta = f"{out_dir}/concat.fasta"
sras_qc_passed_file = f'{out_dir}/sras_qc_passed.txt'
aligned_fasta = f'{out_dir}/ref_aligned.fasta'
aligned_ref_removed_fasta = f'{out_dir}/ref_aligned_samples_only.fasta'

Path(out_dir).mkdir(exist_ok=True)

# read in sras list
target_sras: set[str] = set()
with open(sras_file, 'r') as f:
    target_sras = {l.strip() for l in f.readlines()}

# find input fasta files
all_fastas = glob.glob(f'{fastas_dir}/*HA_cns.fa')

fasta_by_sra = dict()
for fasta in all_fastas:
    fasta_sra = os.path.basename(fasta).replace('_HA_cns.fa', '')
    if fasta_sra in target_sras:
        fasta_by_sra[fasta_sra] = fasta

sras_missing_fastas = target_sras - set(fasta_by_sra.keys())
if len(sras_missing_fastas) > 0:
    print(f'Warning: Fasta files not found for {len(sras_missing_fastas)} SRAs: {sras_missing_fastas}')

seq_by_sra = dict()
for sra, fasta in fasta_by_sra.items():
    with open(fasta, 'r') as f:
        contents = [s for s in SeqIO.parse(f, 'fasta')]
        assert len(contents) == 1 # check that file contains only one sequence
        assert contents[0].id.find(sra) >= 0 # check that file name and contents match
        contents[0].name = sra
        contents[0].id = sra
        contents[0].description = ''
        seq_by_sra[sra] = contents[0]

# QC
# Ambiguous bases
ambiguous_bases = {'N', 'R', 'Y', 'S','W', 'K', 'M', 'B', 'D', 'H', 'V'}
failed_sras = set()
for sra, seq in seq_by_sra.items():
    count_ambig = sum([seq.upper().count(b) for b in ambiguous_bases])
    ambig_frac = count_ambig / len(seq)
    
    if ambig_frac > ambiguous_bases_frac_limit:
        failed_sras.add(sra)

failed_seqs = [seq_by_sra.pop(sra) for sra in failed_sras]
if len(failed_sras) > 0:
    print(f'Warning: {len(failed_sras)} sequences discarded (fraction of ambiguous bases > {ambiguous_bases_frac_limit})')

with open(sras_qc_passed_file, 'w+') as f:
    f.writelines([sra + '\n' for sra in seq_by_sra.keys()])

# concatenate input sequences
with open(concatenated_fasta, 'w+') as f:
    SeqIO.write([seq for seq in seq_by_sra.values()], f, 'fasta')

# Align
mafft_command = [
    'mafft',
    '--auto',
    '--keeplength',
    '--addfragments',
    concatenated_fasta,
    ref_fasta
]
with open(aligned_fasta, 'w+') as f:
    subprocess.run(mafft_command, stdout=f)

# write a version of the MSA with the ref removed
with open(ref_fasta, 'r') as f:
    contents = [s for s in SeqIO.parse(f, 'fasta')]
    assert len(contents) == 1
    ref_seq = contents[0]

with open(aligned_fasta, 'r') as f:
    contents = [s for s in SeqIO.parse(f, 'fasta')]
    aligned_ref_removed = [s for s in contents if s.id != ref_seq.id]

with open(aligned_ref_removed_fasta, 'w+') as f:
    SeqIO.write(aligned_ref_removed, f, 'fasta')
