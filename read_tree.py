from Bio.Nexus.Nexus import Nexus
from Bio.Nexus.Trees import Tree
from Bio.Nexus.Nodes import Node
from Bio import SeqIO
from Bio import Align
from collections import namedtuple
import toytree
import toyplot


Mut = namedtuple('Mut', ['ref', 'pos', 'alt'])

class AaTree:
    def __init__(self, nexus: Nexus, aa_seqs_filename: str):
        assert len(nexus.trees) == 1
        ntree = nexus.trees[0]
        

        # int -> str
        self.name_by_id = dict()
        # int -> int
        self.parent_by_id = dict()
        # int -> {int}
        self.children_by_id = dict()
        # int -> float
        self.parent_dist_by_id = dict()
        self.leaves = set()
        # (ref id: int, other id: int) -> {Mut} 
        self.aa_variants = dict()

        self._build(ntree)

        # str -> int
        self.id_by_name = {v: k for k, v in self.name_by_id.items()}
        # int -> sequence
        self.aa_seq_by_id = dict()

        with open(aa_seqs_filename, 'r') as f:
            aa_seqs = [s for s in SeqIO.parse(f, 'fasta')]
            for s in aa_seqs:
                try:
                    id_ = self.id_by_name[s.name]
                    seq = s.seq
                    self.aa_seq_by_id[id_] = seq
                except KeyError:
                    pass
        assert self.aa_seq_by_id.keys() == self.name_by_id.keys()


    def _build(self, ntree: Tree, nextid: int = 0):
        node = ntree.node(nextid)
        name = node.data.taxon
        if name is None:
            name = node.data.support
        assert name is not None

        self.name_by_id[node.id] = name
        self.parent_by_id[node.id] = node.prev
        self.parent_dist_by_id[node.id] = node.data.branchlength
        children = set(node.get_succ())
        self.children_by_id[node.id] = children

        if len(children) == 0:
            self.leaves.add(node.id)

        for child_id in node.get_succ():
            self._build(ntree, child_id)

    def compare(self, refid: int, qid: int) -> set(Mut):
        try:
            return self.aa_variants[(refid, qid)]
        except KeyError:
            refseq = self.aa_seq_by_id[refid]
            qseq = self.aa_seq_by_id[qid]
            assert len(refseq) == len(qseq)
            muts = set()
            for i in range(len(qseq)):
                    if qseq[i] != refseq[i]:
                        muts.add(Mut(refseq[i], i+1, qseq[i]))
            
            self.aa_variants[(refid, qid)] = muts


    def compare_root_and_leaves(self):
        for leaf_id in self.leaves:
            self.compare(0, leaf_id)

    def compare_parents_and_children(self):
        for cid, pid in self.parent_by_id.items():
            if pid is not None:
                self.compare(pid, cid)


    def to_toytree(self, root_id: int = 0) -> toytree.Node:
        n = toytree.Node(
                self.name_by_id[root_id],
                dist=self.parent_dist_by_id[root_id]
            )
        n

        for cid in self.children_by_id[root_id]:
            n._add_child(self.to_toytree(cid))
        return n

    def get_new_muts_by_name(self, name: str) -> set(Mut):
        id_ = self.id_by_name[name]
        parent_id = self.parent_by_id[id_]
        try:
            return self.aa_variants[(parent_id, id_)]
        except KeyError:
            return set()




with open('out/ancestral/annotated_tree.nexus', 'r') as f:
    direct_toytree = toytree.tree(f.read())

with open('out/ancestral/annotated_tree.nexus', 'r') as f:
    nexus = Nexus()
    nexus.read(f)


# todo: this should come from a real source.
seq_file = 'out/playground/amino_acids.fa'
aatree = AaTree(nexus, seq_file)

aatree.compare_root_and_leaves()
aatree.compare_parents_and_children()

new_aa_muts_by_name = {name: {str(m) for m in aatree.get_new_muts_by_name(name)} for name in aatree.id_by_name.keys()}
count_new_aa_muts_by_name = {k: len(v) for k, v in new_aa_muts_by_name.items()}
print(count_new_aa_muts_by_name)

direct_toytree = direct_toytree.set_node_data(feature='count_new_aa_muts', data=count_new_aa_muts_by_name, default=0)
print(direct_toytree.get_nodes('SRR29182385'))
print(direct_toytree)
print(direct_toytree.features)

print(direct_toytree.get_node_data('count_new_aa_muts'))

# we got the names for free! 
for n in direct_toytree.traverse():
    # print(n)
    pass

# toy_root = aatree.to_toytree()
# tree = toytree.tree(toy_root)

canvas, axes, mark = direct_toytree.draw(tip_labels=True, height=8000, width=4000, tip_labels_style={ 
        "font-size": "8px", 
    });
toytree.save(canvas, "/tmp/drawing.svg")

