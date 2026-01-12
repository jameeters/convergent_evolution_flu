from Bio import SeqIO
from collections import namedtuple
import toytree
import toyplot
from typing import NamedTuple
from pprint import pprint


class Mut(NamedTuple):
    ref: str
    pos: int
    alt: str

    def __str__(self):
        return f'{self.ref}{self.pos}{self.alt}'

    def __repr__(self):
        return f'Mut({str(self)})'


def compare_aa_seqs(base: toytree.Node, sample: toytree.Node) -> set(Mut):
    base_seq = base.aa_seq
    sample_seq = sample.aa_seq
    assert len(base_seq) == len(sample_seq)
    muts = set()
    for i in range(len(sample_seq)):
            if sample_seq[i] != base_seq[i] and sample_seq[i] != 'X':
                muts.add(Mut(base_seq[i], i+1, sample_seq[i]))
    return muts


def compare_nt_seq(base: toytree.Node, sample: toytree.Node) -> set(Mut):
    base_seq = base.nt_seq
    sample_seq = sample.nt_seq
    assert len(base_seq) == len(sample_seq)
    muts = set()
    for i in range(len(sample_seq)):
            if sample_seq[i] != base_seq[i] and sample_seq[i] != 'N':
                muts.add(Mut(base_seq[i], i+1, sample_seq[i]))
    return muts


with open('out/ancestral/annotated_tree.nexus', 'r') as f:
    direct_toytree = toytree.tree(f.read())

nt_seqs_filename = 'out/translation/coding_seqs.fasta'
nt_seq_by_name = dict()
with open(nt_seqs_filename, 'r') as f:
            nt_seqs = [s for s in SeqIO.parse(f, 'fasta')]
            for s in nt_seqs:
                seq = s.seq
                nt_seq_by_name[s.name] = seq

direct_toytree = direct_toytree.set_node_data(feature='nt_seq', data=nt_seq_by_name)

# todo: this should come from a real source.
aa_seqs_filename = 'out/playground/amino_acids.fasta'

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
new_nt_muts_by_name = dict()
for node in direct_toytree:
    parent = node.up
    aa_muts = set()
    nt_muts = set()
    if parent is not None:
        aa_muts = compare_aa_seqs(base=parent, sample=node)
        nt_muts = compare_nt_seq(base=parent, sample=node)
    new_aa_muts_by_name[node.name] = aa_muts
    new_nt_muts_by_name[node.name] = nt_muts



count_new_aa_muts_by_name = {k: len(v) for k, v in new_aa_muts_by_name.items()}
direct_toytree = direct_toytree.set_node_data(feature='count_new_aa_muts', data=count_new_aa_muts_by_name, default=0)
direct_toytree = direct_toytree.set_node_data(feature='new_aa_muts', data=new_aa_muts_by_name, default=set())
direct_toytree = direct_toytree.set_node_data(feature='new_nt_muts', data=new_nt_muts_by_name, default=set())

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
            putative = putative.up
        try:
            leaf_aa_muts_from_root_emergence[mut][leaf.name] = putative.name
        except KeyError:
            leaf_aa_muts_from_root_emergence[mut] = {leaf.name: putative.name}

distinct_emergences = {mut: len(set(v.values())) for mut, v in leaf_aa_muts_from_root_emergence.items()}

# for mut, count in distinct_emergences.items():
#     if count > 1:
#         print(mut)
#         pprint(leaf_aa_muts_from_root_emergence[mut])

new_a172t = {node.name: True for node in direct_toytree if Mut('A', 172, 'T') in node.new_aa_muts}
direct_toytree = direct_toytree.set_node_data(feature='has_A172T', data=new_a172t, inherit=True, default=False)

# canvas, axes, mark = direct_toytree.draw(
#     tip_labels=True,
#     height=8000,
#     width=4000,
#     tip_labels_style={ "font-size": "8px",},
#     node_colors=('count_new_aa_muts', 'BlueRed'),
#     node_sizes=10,
#     node_mask=False,
#     node_labels='count_new_aa_muts',
#     node_labels_style={'font-size': '8px', 'fill': '#66aa66'},
#     node_style={'stroke': None},
# )

canvas, axes, mark = direct_toytree.draw(
    tip_labels=True,
    height=8000,
    width=4000,
    tip_labels_style={ "font-size": "8px",},
    node_colors=('has_A172T', 'BlueRed'),
    node_sizes=5,
    node_mask=False,
    node_labels=None,
    node_labels_style={'font-size': '8px', 'fill': '#66aa66'},
    node_style={'stroke': None},
)

toytree.save(canvas, "/tmp/drawing.svg")



print(direct_toytree.features)
for node in direct_toytree:
    try:
        if len(node.new_nt_muts) > 0:
            print(f'{node.name} :: {node.new_nt_muts}')
    except AttributeError:
        pass