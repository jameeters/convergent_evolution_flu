from Bio import SeqIO
from collections import namedtuple
import toytree
import toyplot
from typing import NamedTuple
from pprint import pprint
from Bio.Seq import Seq
from utils.translation import translate_seq_codonwise
import pandas as pd


# Inputs
tree_file = 'out/ancestral/annotated_tree.nexus'
nt_seqs_filename = 'out/translation/coding_seqs.fasta'
aa_seqs_filename = 'out/translation/amino_acids.fasta'

# Outputs
codon_mutation_counts_output_tsv = 'out/codon_mutation_counts_bovine_ha.tsv'
dn_ds_output_tsv = 'out/dn_ds_by_site_bovine_ha.tsv'

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

# ------------------------------------------------
# count syn and non-syn mutations per site
# ------------------------------------------------
mutation_counts = dict()
for mut_set in new_codon_muts_by_name.values():
    for mut in mut_set:
        
        # Lists are not hashable, so this will force an error if values remain unset. 
        incr_count_name = ['unset']
        other_count_name = ['unset']
        if mut.is_synonymous():
            incr_count_name = 'count_s'
            other_count_name = 'count_n'
        else:
            incr_count_name = 'count_n'
            other_count_name = 'count_s'

        try:
            mutation_counts[mut.pos_aa][incr_count_name] += 1
        except KeyError:
            mutation_counts[mut.pos_aa] = {'pos_aa': mut.pos_aa, incr_count_name: 1, other_count_name: 0}

mutation_counts_df = pd.DataFrame(data=mutation_counts.values())

# copy Praneeth's formula to avoid division by 0
mutation_counts_df['ratio'] = mutation_counts_df['count_n'] / (mutation_counts_df['count_s'] + 1)

with open(dn_ds_output_tsv, 'w+') as f:
    mutation_counts_df.to_csv(f, sep='\t', index=False)

# ------------------------------------------------
# Emergence of aa mutations in leaves wrt root. 
# output from this is not currently used
# ------------------------------------------------
leaf_aa_muts_from_root_by_name = dict()
for leaf in direct_toytree.treenode.get_leaves():
    leaf_aa_muts_from_root_by_name[leaf.name] = compare_aa_seqs(base=direct_toytree.treenode, sample=leaf)
direct_toytree = direct_toytree.set_node_data(feature='aa_muts_wrt_root', data=leaf_aa_muts_from_root_by_name, default=None)

# mut -> leaf name -> name of node where mut emerged
# {Mut: {str: str}}
leaf_aa_muts_from_root_emergence = dict()

for leaf in direct_toytree.treenode.get_leaves():
    for mut in leaf.aa_muts_wrt_root:
        putative = leaf
        while mut not in putative.new_aa_muts:

            # hacky fix to account for multiple mutations at the same site.
            # EG: at leaf, mut wrt root is K228I and mut wrt parent is R228I. 
            if len([pm for pm in putative.new_aa_muts if pm.pos == mut.pos and pm.alt == mut.alt]) > 0:
                break

            putative = putative.up

            if putative == direct_toytree.treenode:
                print(f'Got all the way to root: mutation {mut} observed in leaf {leaf.name}')
                break
        try:
            leaf_aa_muts_from_root_emergence[mut][leaf.name] = putative.name
        except KeyError:
            leaf_aa_muts_from_root_emergence[mut] = {leaf.name: putative.name}

distinct_emergences = {mut: len(set(v.values())) for mut, v in leaf_aa_muts_from_root_emergence.items()}


canvas, axes, mark = direct_toytree.draw(
    tip_labels=True,
    height=8000,
    width=4000,
    tip_labels_style={ "font-size": "8px",},
    node_colors=('count_new_aa_muts', 'BlueRed'),
    node_sizes=10,
    node_mask=False,
    node_labels='count_new_aa_muts',
    node_labels_style={'font-size': '8px', 'fill': '#66aa66'},
    node_style={'stroke': None},
)

toytree.save(canvas, "/tmp/drawing.svg")
