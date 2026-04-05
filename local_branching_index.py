from Bio import SeqIO, AlignIO
from collections import namedtuple
import toytree
import pandas as pd
import argparse
import numpy as np
import baltic as bt
import matplotlib as mpl
from matplotlib import pyplot as plt



def main(args):

    tree = read_tree(args.treefile)

    # todo: it's possible that this should be using the divergence tree
    # or some other pure genetic tree rather than the time tree.
    # I have confirmed that the timetree.nexus and divergence_tree.nexus have 
    # exactly the same topological structure and node names:
    #
    # cat out_avian_b3.13_ha/treetime/divergence_tree.nexus | grep -E '^ Tree.*' | cut -d '=' -f 2- | grep -oE '([()]|(SRR|NODE_)\d+)' | tr -d '\n' > /tmp/foo
    # cat out_avian_b3.13_ha/treetime/timetree.nexus | grep -E '^ Tree.*' | cut -d '=' -f 2- | grep -oE '([()]|(SRR|NODE_)\d+)' | tr -d '\n' > /tmp/bar
    # md5sum /tmp/{foo,bar}

    # calculate tau (magic number from neher et al 2014)
    tau = 0.0625 

    # Neher's impl uses time_scale = pairwise dist / 2, and also multiplies tau by 2
    tau *= 2
    time_scale = calc_mean_pairwise_dist(args.sample_msa_fasta) / 2

    # value from paper, but not using this yet. shim_node_dist = 1e-7 
    # I'm not sure what problems occur if branch len = 0, I don'e think it's 
    # a denominator anywhere...
    shim_node_dist = 1e-7


    calc_lbi(tree, tau, time_scale, shim_node_dist)

    #plot_tree(tree)

    with open(f'{args.outdir}/local_branching_index_shim.tsv', 'w+') as f:
        lines = ['node\tlbi\n']
        lines += [f'{n.name}\t{n.lbi}\n' for n in tree]
        f.writelines(lines)


def calc_mean_pairwise_dist(msa_fasta: str):
    aln = AlignIO.read(msa_fasta, 'fasta')
    n_nt = len(aln) * aln.get_alignment_length()
    nt_freqs = np.array([v for v in aln.alignment.frequencies.values()]) / n_nt
    
    return np.sum(nt_freqs * (1 - nt_freqs), axis=0).mean()



def calc_lbi(tree, tau, time_scale, shim_node_dist):

    # Traverse up from tips to root, calculating up-messages
    for node in tree.traverse(strategy='idxorder'):
        m_up = 0
        bl = (node.dist + shim_node_dist) / time_scale
        for child in node.children:
            m_up += child.m_up
        m_up *= np.exp(-bl / tau)
        m_up += tau * (1 - np.exp(-bl / tau))
        node.m_up = m_up


    # traverse from root to tips, calculating down-messages
    tree.treenode.m_down = 0
    for node in tree.traverse(strategy='levelorder'):
        child_up_messages_sum = sum([child.m_up for child in node.children])
        for child in node.children:
            bl = (child.dist + shim_node_dist) / time_scale
            # exclude this child from the sum
            m_down = child_up_messages_sum - child.m_up 
            m_down += node.m_down
            m_down *= np.exp(-bl / tau)
            m_down += tau * (1 - np.exp(-bl / tau))
            child.m_down = m_down

    # Traverse the tree one more time to calculate all the final LBI values
    for node in tree.traverse(strategy='levelorder'):
        node.lbi = node.m_down + sum([c.m_up for c in node.children])


def plot_tree(tree):
    with open('/tmp/tmp_tree_mcgee.newick', 'w+') as f:
        tree.write('/tmp/tmp_tree_mcgee.newick', features=['m_up', 'm_down', 'lbi'])
        bttree = bt.loadNewick(f)

    min_lbi = 200
    max_lbi = 0
    for node in tree.traverse(strategy='levelorder'):
        min_lbi = min([min_lbi, node.lbi])
        max_lbi = max([max_lbi, node.lbi])

    cmap = mpl.cm.viridis

    # x_attr = lambda k: k.traits['date']
    c_func = lambda k: cmap(k.traits['lbi']/max_lbi)

    fig, ax = plt.subplots(figsize=(10,20),facecolor='w')

    norm = mpl.colors.Normalize(vmin=min_lbi, vmax=max_lbi)

    bttree.plotTree(ax, colour='k')
    bttree.plotPoints(ax, colour=c_func)
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical', label='LBI')
    ax.plot()

    plt.show()

def read_tree(tree_file: str):
    with open(tree_file, 'r') as f:
        return toytree.tree(f.read())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('treefile')
    parser.add_argument('sample_msa_fasta')
    parser.add_argument('outdir')
    
    args = parser.parse_args()
    main(args)
