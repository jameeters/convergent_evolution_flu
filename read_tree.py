from Bio.Nexus.Nexus import Nexus
from Bio.Nexus.Trees import Tree
from Bio.Nexus.Nodes import Node
from Bio import SeqIO


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

        self._build(ntree)
        self.id_by_name = {v: k for k, v in self.name_by_id.items()}

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
        self.children_by_id[node.id] = set(node.get_succ())

        for child_id in node.get_succ():
            self._build(ntree, child_id)


def add_to_map(m: dict, id: int, tree: Tree) -> dict:
    node = tree.node(id)
    name = node.data.taxon
    if name is None:
        name = node.data.support
    assert name is not None

    m[node.id] = name

    for child_id in node.get_succ():
        m = add_to_map(m, child_id, tree)
    return m

def map_node_ids(tree: Tree) -> dict:
    id_map = dict()
    id_map = add_to_map(id_map, 0, tree)
    return id_map


with open('out/ancestral/annotated_tree.nexus', 'r') as f:
    nexus = Nexus()
    nexus.read(f)



seq_file = 'out/playground/amino_acids.fa'
aatree = AaTree(nexus, seq_file)

print(aatree.aa_seq_by_id)



