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

# Call sample and ancestral mutations with respect to DMS reference

def main(args):
    aa_seqs = read_fasta(args.aa_seqs_fasta)
    ref_seq = get_ref_aa_seq(args.ref_fasta)

    mutations = call_mutations(ref_seq, aa_seqs)
    
    out_filename = f'{args.out_dir}/mutations_wrt_ref.tsv'
    with open(out_filename, 'w+') as f:
        mutations.to_csv(f, sep='\t', index=False)

def call_mutations(ref_seq, aa_seqs) -> pd.DataFrame:
    mutations = []
    for name, samp_seq in aa_seqs.items():
        assert len(samp_seq) == len(ref_seq)
        for i, samp_aa in enumerate(samp_seq):
            ref_aa = ref_seq[i]
            if samp_aa != ref_aa and samp_aa != 'X':
                mutations.append({
                    'node': name,
                    'pos_aa': i + 1,
                    'ref_aa': ref_aa,
                    'alt_aa': samp_aa
                    })

    mutations_df = pd.DataFrame(data=mutations)
    return mutations_df


def read_fasta(filename: str) -> Dict[str, Seq]:
    aa_seq_by_name = dict()
    with open(filename, 'r') as f:
        aa_seqs = [s for s in SeqIO.parse(f, 'fasta')]
        for s in aa_seqs:
            seq = s.seq
            aa_seq_by_name[s.name] = seq
    return aa_seq_by_name

def get_ref_aa_seq(filename: str) -> Seq:
    ref_seq_dict = read_fasta(filename)
    assert len(ref_seq_dict) == 1
    nt_seq = [v for v in ref_seq_dict.values()][0]
    return nt_seq.translate()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # parse args
    parser.add_argument('ref_fasta')
    parser.add_argument('aa_seqs_fasta')

    parser.add_argument('out_dir')

    args = parser.parse_args()
    main(args)