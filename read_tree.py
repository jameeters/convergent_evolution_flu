from Bio import SeqIO
from collections import namedtuple
import toytree
import toyplot
from typing import NamedTuple
from pprint import pprint
from Bio.Seq import Seq
from utils.translation import translate_seq_codonwise
import pandas as pd
import argparse


parser = argparse.ArgumentParser()

# Inputs
parser.add_argument('translation_out_dir')
parser.add_argument('treetime_out_dir')


# Outputs
parser.add_argument('out_dir')

args = parser.parse_args()

# Inputs
tree_file = f'{args.treetime_out_dir}/timetree.nexus'
nt_seqs_filename = f'{args.translation_out_dir}/coding_seqs.fasta'
aa_seqs_filename = f'{args.translation_out_dir}/amino_acids.fasta'
node_dates_filename = f'{args.treetime_out_dir}/dates.tsv'

# Outputs
codon_mutation_counts_output_tsv = f'{args.out_dir}/codon_mutation_counts.tsv'
dated_mutations_tsv = f'{args.out_dir}/dated_mutations.tsv'

class CodonMut(NamedTuple):
    ref_nt: str
    alt_nt: str
    pos_start_nt: int
    ref_aa: str
    alt_aa: str
    pos_aa: int

    def is_synonymous(self) -> bool:
        return self.ref_aa == self.alt_aa

    def to_dict(self, count: int = -1) -> dict:
        return {
            'ref_nt': self.ref_nt,
            'alt_nt': self.alt_nt,
            'pos_start_nt': self.pos_start_nt,
            'ref_aa': self.ref_aa,
            'alt_aa': self.alt_aa,
            'pos_aa': self.pos_aa,
            'count': count
        }

    def __str__(self):
        return f'{self.ref_nt}({self.ref_aa}) {self.pos_start_nt}({self.pos_aa}) {self.alt_nt}({self.alt_aa})'

    def __repr__(self):
        return f'CodonMut({str(self)})'

class Mut(NamedTuple):
    ref: str
    pos: int
    alt: str

    def __str__(self):
        return f'{self.ref}{self.pos}{self.alt}'

    def __repr__(self):
        return f'Mut({str(self)})'


Count_n_s = namedtuple('Count_n_s', ['count_n', 'count_s', 'ratio'])


def compare_aa_seqs(base: toytree.Node, sample: toytree.Node) -> set(Mut):
    base_seq = base.aa_seq
    sample_seq = sample.aa_seq
    assert len(base_seq) == len(sample_seq)
    muts = set()
    for i in range(len(sample_seq)):
            if sample_seq[i] != base_seq[i] and sample_seq[i] != 'X':
                muts.add(Mut(base_seq[i], i+1, sample_seq[i]))
    return muts


def compare_codons(base: toytree.Node, sample: toytree.Node) -> set(Mut):
    base_seq = base.nt_seq
    sample_seq = sample.nt_seq
    assert len(base_seq) == len(sample_seq)
    muts = set()

    for i in range(0, len(sample_seq), 3):
        ref = base_seq[i:i+3]
        alt = sample_seq[i:i+3]

        if ref != alt and alt != 'nnn':
            ref_aa = translate_seq_codonwise(Seq(ref))
            alt_aa = translate_seq_codonwise(Seq(alt))
            if alt_aa == 'X':
                continue
            pos_aa = (i/3) + 1
            assert pos_aa == int(pos_aa)
            pos_aa = int(pos_aa)
            muts.add(CodonMut(ref, alt, i+1, ref_aa, alt_aa, pos_aa))
    return muts

def check_codons_and_amino_acids(node: toytree.Node):
    for cm in node.new_codon_muts:
        if not cm.is_synonymous():
            assert Mut(ref=cm.ref_aa, alt=cm.alt_aa, pos=cm.pos_aa) in node.new_aa_muts


def read_dates():
    with open(node_dates_filename, 'r') as f:
        dates = pd.read_csv(f, sep='\t')
    return dates



with open(tree_file, 'r') as f:
    direct_toytree = toytree.tree(f.read())


nt_seq_by_name = dict()
with open(nt_seqs_filename, 'r') as f:
            nt_seqs = [s for s in SeqIO.parse(f, 'fasta')]
            for s in nt_seqs:
                if s.name == 'HA|PP755589.1|A/cattle/Texas/24-008749-003/2024(H5N1)':
                    continue
                seq = s.seq
                nt_seq_by_name[s.name] = seq

direct_toytree = direct_toytree.set_node_data(feature='nt_seq', data=nt_seq_by_name)


aa_seq_by_name = dict()
with open(aa_seqs_filename, 'r') as f:
            aa_seqs = [s for s in SeqIO.parse(f, 'fasta')]
            for s in aa_seqs:
                if s.name == 'HA|PP755589.1|A/cattle/Texas/24-008749-003/2024(H5N1)':
                    continue
                seq = s.seq
                aa_seq_by_name[s.name] = seq


direct_toytree = direct_toytree.set_node_data(feature='aa_seq', data=aa_seq_by_name)

# compare parents to children
new_aa_muts_by_name = dict()
new_codon_muts_by_name = dict()
for node in direct_toytree:
    parent = node.up
    aa_muts = set()
    codon_muts = set()
    if parent is not None:
        aa_muts = compare_aa_seqs(base=parent, sample=node)
        codon_muts = compare_codons(base=parent, sample=node)
    new_aa_muts_by_name[node.name] = aa_muts
    new_codon_muts_by_name[node.name] = codon_muts


count_new_aa_muts_by_name = {k: len(v) for k, v in new_aa_muts_by_name.items()}
direct_toytree = direct_toytree.set_node_data(feature='count_new_aa_muts', data=count_new_aa_muts_by_name, default=0)
direct_toytree = direct_toytree.set_node_data(feature='new_aa_muts', data=new_aa_muts_by_name, default=set())
direct_toytree = direct_toytree.set_node_data(feature='new_codon_muts', data=new_codon_muts_by_name, default=set())

# spot check that two different translation methods match
for node in direct_toytree:
    check_codons_and_amino_acids(node)

# ---------------------------------
# output full codon mutation data
# ---------------------------------
codon_muts_counts = dict()

for mut_set in new_codon_muts_by_name.values():
    for mut in mut_set:
        try:
            codon_muts_counts[mut] += 1
        except KeyError:
            codon_muts_counts[mut] = 1
temp = [mut.to_dict(count) for mut, count in codon_muts_counts.items()]
codon_muts_counts_df = pd.DataFrame(data=temp)


with open(codon_mutation_counts_output_tsv, 'w+') as f:
    codon_muts_counts_df.to_csv(f, sep='\t', index=False)

# dated mutation data
dates = read_dates()
dates = dates.rename(columns={'#node': 'node'})

temp = {node: [mut.to_dict() for mut in mutset] for node, mutset in new_codon_muts_by_name.items()}
for node, mutlist in temp.items():
    for mut in mutlist:
        mut['node'] = node

temp2 = []
for mutlist in temp.values():
    temp2 += mutlist

dated_muts_df = pd.DataFrame(data=temp2).merge(dates, how='left', on = 'node', validate="m:1")

with open(dated_mutations_tsv, 'w+') as f:
    dated_muts_df.to_csv(f, sep='\t', index=False)
