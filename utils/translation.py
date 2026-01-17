from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data.CodonTable import TranslationError


def translate_seq_codonwise(seq: Seq) -> Seq:
    amino_acids = Seq('')
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        try:
            amino_acids += codon.translate()
        except TranslationError:
            if '-' in codon:
                amino_acids += 'X'
    return amino_acids

def translate_codonwise(seqrec: SeqRecord) -> SeqRecord:
    amino_acids = translate_seq_codonwise(seqrec.seq)
    return SeqRecord(amino_acids, id=seqrec.id, name=seqrec.name, description=seqrec.description)