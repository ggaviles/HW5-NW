# Importing Dependencies
import pytest
from Bio import Align
from Bio.Align import substitution_matrices
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    needle_wunsch = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10.0, -1.0)
    alignment_score, seq1_align, seq2_align = needle_wunsch.align(seq1, seq2)

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    # Had to change this to -11 instead of -10 because Align.PairwiseAligner() does not add gap_open and gap_extend
    # for one initial gap, while my model does. To compensate for that difference, I changed the open_gap_score to be
    # my open gap score + extend gap score (-10 + -1), which is -11. Since there is a singular gap in the alignment
    # between seq1 and seq2, I think this is a reasonable change.
    aligner.open_gap_score = -11

    # For some reason, I get a score of
    for alignment in aligner.align(seq1, seq2):
        assert alignment_score == alignment.score

    # Alignment is correct
    assert seq1_align == "MYQR"
    assert seq2_align == "M-QR"

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    needle_wunsch = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10.0, -1.0)
    alignment_score, seq3_align, seq4_align = needle_wunsch.align(seq3, seq4)

    assert alignment_score == 17
    assert seq3_align == 'MAVHQLIRRP'
    assert seq4_align == 'M---QLIRHP'




