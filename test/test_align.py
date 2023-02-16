# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    needle_wunsch = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10.0, -1.0)


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




