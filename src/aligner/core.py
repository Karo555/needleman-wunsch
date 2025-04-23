from typing import List, Tuple
from src.aligner.models import Sequence


def build_score_matrix(
    seq1: Sequence, seq2: Sequence, match: int, mismatch: int, gap: int
) -> List[List[int]]:
    """
    Build and return the scoring matrix for global alignment using
    the Needleman–Wunsch algorithm.

    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for a match
    :param mismatch: Score for a mismatch
    :param gap: Gap penalty (negative)
    :return: A (len(seq1)+1) x (len(seq2)+1) matrix of scores
    """
    n = len(seq1)
    m = len(seq2)
    # initialize matrix with zeros
    matrix: List[List[int]] = [[0] * (m + 1) for _ in range(n + 1)]

    # initialize first column with gap penalties
    for i in range(1, n + 1):
        matrix[i][0] = matrix[i - 1][0] + gap

    # initialize first row with gap penalties
    for j in range(1, m + 1):
        matrix[0][j] = matrix[0][j - 1] + gap

    # fill the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            char1 = seq1.sequence[i - 1]
            char2 = seq2.sequence[j - 1]
            # match or mismatch
            if char1 == char2:
                diag = matrix[i - 1][j - 1] + match
            else:
                diag = matrix[i - 1][j - 1] + mismatch
            # gap in seq2 (up)
            up = matrix[i - 1][j] + gap
            # gap in seq1 (left)
            left = matrix[i][j - 1] + gap
            # choose max
            matrix[i][j] = max(diag, up, left)

    return matrix


def traceback(
    matrix: List[List[int]],
    seq1: Sequence,
    seq2: Sequence,
    match: int,
    mismatch: int,
    gap: int,
) -> Tuple[str, str]:
    """
    Perform a traceback through the scoring matrix to recover one optimal alignment.

    :param matrix: Scoring matrix from build_score_matrix
    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for a match
    :param mismatch: Score for a mismatch
    :param gap: Gap penalty
    :return: A tuple of aligned strings (with '-' for gaps)
    """
    i, j = len(seq1), len(seq2)
    aligned1: List[str] = []
    aligned2: List[str] = []

    # Traceback from bottom-right to top-left
    while i > 0 or j > 0:
        # when both sequences have characters remaining
        if i > 0 and j > 0:
            char1 = seq1.sequence[i - 1]
            char2 = seq2.sequence[j - 1]
            # score if coming from diagonal
            if char1 == char2:
                score_diag = matrix[i - 1][j - 1] + match
            else:
                score_diag = matrix[i - 1][j - 1] + mismatch
            # check if diagonal is optimal
            if matrix[i][j] == score_diag:
                aligned1.append(char1)
                aligned2.append(char2)
                i -= 1
                j -= 1
                continue
        # check if coming from up (gap in seq2)
        if i > 0 and matrix[i][j] == matrix[i - 1][j] + gap:
            aligned1.append(seq1.sequence[i - 1])
            aligned2.append("-")
            i -= 1
            continue
        # otherwise, must be coming from left (gap in seq1)
        if j > 0 and matrix[i][j] == matrix[i][j - 1] + gap:
            aligned1.append("-")
            aligned2.append(seq2.sequence[j - 1])
            j -= 1
            continue
        # safety break (should not happen)
        break

    # Reverse to get the correct order
    aligned1.reverse()
    aligned2.reverse()
    return "".join(aligned1), "".join(aligned2)
