# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment

        gap_open_penalty = self.gap_open
        gap_ext_penalty = self.gap_extend

        sub_dict = self.sub_dict
        # create matrices for alignment scores, gaps, and backtracing

        # Determine length of each sequence and store in n or m
        n = len(seqA)
        m = len(seqB)

        # Generate matrix of zeros to store scores
        self._align_matrix = np.zeros((m + 1, n + 1))
        self._gapA_matrix = np.zeros((m + 1, n + 1))
        self._gapB_matrix = np.zeros((m + 1, n + 1))

        # TODO: Implement global alignment here

        # Initialize empty matrices
        score_mat = np.zeros((m + 1, n + 1))

        # Initialize matrix values

        for i in range(0, m + 1):
            score_mat[i][0] = gap_open_penalty * i
            self._align_matrix[i][0] = gap_open_penalty * i
            self._gapA_matrix[i][0] = np.NINF
            self._gapB_matrix[i][0] = -(gap_open_penalty + gap_ext_penalty * i)

        for j in range(0, n + 1):
            score_mat[0][j] = gap_open_penalty * j
            self.align_matrix[0][j] = gap_open_penalty * j
            self._gapA_matrix[0][j] = -(gap_open_penalty + gap_ext_penalty * j)
            self._gapB_matrix[0][j] = np.NINF

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                score_mat[i][j] = (self._align_matrix[i - 1][j - 1] + sub_dict[(seqA[j - 1], seqB[i - 1])]) + max(
                    score_mat[i - 1][j - 1], self._gapA_matrix[i - 1][j - 1], self._gapB_matrix[i - 1][j - 1])
                self._gapA_matrix[i][j] = max((score_mat[i][j - 1] + (gap_open_penalty + gap_ext_penalty)),
                                     (self._gapA_matrix[i][j - 1] + gap_open_penalty),
                                     ((self._gapB_matrix[i][j - 1]) + (gap_open_penalty + gap_ext_penalty)))
                self._gapB_matrix[i][j] = max((score_mat[i - 1][j] + (gap_open_penalty + gap_ext_penalty)),
                                     (self._gapA_matrix[i - 1][j] + (gap_open_penalty + gap_ext_penalty)),
                                     (self._gapB_matrix[i - 1][j] + gap_open_penalty))

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        seqA = self._seqA
        seqB = self._seqB
        sub_dict = self.sub_dict
        alignA = self.seqA_align
        alignB = self.seqB_align

        gap_penalty = self.gap_open

        i = len(seqB)
        j = len(seqA)

        while i > 0 and j > 0:  # end touching the top or the left edge
            score_current = self._align_matrix[i][j]
            score_diagonal = self._align_matrix[i - 1][j - 1]
            score_up = self._align_matrix[i][j - 1]
            score_left = self._align_matrix[i - 1][j]

            # Check to figure out which cell the current score was calculated from,
            # then update i and j to correspond to that cell.
            if score_current == score_diagonal + sub_dict[(seqA[j - 1], seqB[i - 1])]:
                alignA += seqA[j - 1]
                alignB += seqB[i - 1]
                i -= 1
                j -= 1
            elif score_current == score_up + gap_penalty:
                alignA += seqA[j - 1]
                alignB += '-'
                j -= 1
            elif score_current == score_left + gap_penalty:
                alignA += '-'
                alignB += seqB[i - 1]
                i -= 1

        # Finish tracing up to the top left cell
        while i > 0:
            alignA += seqA[i - 1]
            alignB += '-'
            j -= 1
        while j > 0:
            alignA += '-'
            alignB += seqB[j - 1]
            i -= 1

        # Traversed the score matrix from the bottom right so the two sequences will be reversed.
        # Reverse the sequences.
        alignA = alignA[::-1]
        alignB = alignB[::-1]

        # Update seq alignment values
        self.seqA_align = alignA
        self.seqB_align = alignB

        self.alignment_score = self._align_matrix.sum()

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
