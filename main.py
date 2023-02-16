# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    needle_wunsch = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10.0, -1.0)

    # Align Gallus gallus to Homo sapiens
    align_gallus_score, hs_gallus_align, gallus_hs_align = needle_wunsch.align(hs_seq, gg_seq)

    # Align Mus musculus to Homo sapiens
    align_mus_score, hs_mus_align, mus_hs_align = needle_wunsch.align(hs_seq, mm_seq)

    # Align Balaeniceps rex to Homo sapiens
    align_bal_score, hs_bal_align, bal_hs_align = needle_wunsch.align(hs_seq, br_seq)

    # Align Tursiops truncatus to Homo sapiens
    align_turs_score, hs_turs_align, turs_hs_align = needle_wunsch.align(hs_seq, tt_seq)

    # TODO print all of the alignment scores between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print('Gallus gallus alignment score: ', align_gallus_score)
    print('Mus musculus alignment score: ', align_mus_score)
    print('Balaeniceps rex alignment score: ', align_bal_score)
    print('Tursiops truncatus alignment score: ', align_turs_score)

if __name__ == "__main__":
    main()
