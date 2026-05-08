from itertools import product

import numpy as np
from Bio.Align import substitution_matrices
from numba import njit, prange
from numba_progress import ProgressBar

from config import Stage1Config
from preprocessing import PreprocessedDataset


def _load_scoring_matrix(scoring_matrix_name: str) -> np.ndarray:
    str_scoring_matrix = substitution_matrices.load(scoring_matrix_name)
    scoring_matrix = np.full((26, 26), -1.0)

    for i, j in product(range(26), range(26)):
        char_i = chr(ord("A") + i)
        char_j = chr(ord("A") + j)
        try:
            scoring_matrix[i, j] = str_scoring_matrix[char_i, char_j]  # pyright: ignore[reportArgumentType,reportCallIssue]
        except IndexError:
            pass

    return scoring_matrix


@njit(nogil=True)
def _smith_waterman(
    seq_a: np.ndarray,
    seq_b: np.ndarray,
    scoring_matrix: np.ndarray,
    gap_open: float,
    gap_extend: float,
) -> float:
    len_a = len(seq_a)
    len_b = len(seq_b)

    prev_row = np.zeros(len_b + 1, dtype=np.float64)
    curr_row = np.zeros(len_b + 1, dtype=np.float64)
    gap_in_b_score = np.full(len_b + 1, float("-inf"), dtype=np.float64)

    best_score = 0.0
    for pos_a in range(1, len_a + 1):
        gap_in_a_score = float("-inf")
        residue_a = seq_a[pos_a - 1]
        for pos_b in range(1, len_b + 1):
            residue_b = seq_b[pos_b - 1]
            gap_in_a_score = max(
                curr_row[pos_b - 1] - gap_open,
                gap_in_a_score - gap_extend,
            )
            gap_in_b_score[pos_b] = max(
                prev_row[pos_b] - gap_open,
                gap_in_b_score[pos_b] - gap_extend,
            )

            diag_score = prev_row[pos_b - 1] + scoring_matrix[residue_a, residue_b]
            cell_score = diag_score

            if gap_in_a_score > cell_score:
                cell_score = gap_in_a_score
            if gap_in_b_score[pos_b] > cell_score:
                cell_score = gap_in_b_score[pos_b]
            if cell_score < 0.0:
                cell_score = 0.0

            curr_row[pos_b] = cell_score
            if cell_score > best_score:
                best_score = cell_score

        prev_row, curr_row = curr_row, prev_row

    return best_score


@njit(parallel=True, nogil=True)
def _jitted_score_protein_pairs(
    proteins: np.ndarray,
    lengths: np.ndarray,
    pairs: np.ndarray,
    scoring_matrix: np.ndarray,
    gap_open: float,
    gap_extend: float,
    progress: ProgressBar,
) -> np.ndarray:
    scores = np.empty(len(pairs), dtype=np.float64)

    for pair_i in prange(len(pairs)):
        pair = pairs[pair_i]
        protein_a_i, protein_b_i = pair
        protein_a, protein_b = (
            proteins[protein_a_i, : lengths[protein_a_i]],
            proteins[protein_b_i, : lengths[protein_b_i]],
        )

        scores[pair_i] = _smith_waterman(
            protein_a, protein_b, scoring_matrix, gap_open, gap_extend
        )

        progress.update(1)

    return scores


def score_protein_pairs(
    pair_list: list[tuple[int, int]],
    preprocessed_dataset: PreprocessedDataset,
    config: Stage1Config,
) -> dict[tuple[int, int], float]:
    scoring_matrix = _load_scoring_matrix(config.alignment_scoring_matrix)

    proteins_list = []
    for protein in preprocessed_dataset.proteins:
        proteins_list.append([ord(char) - ord("A") for char in protein.protein_str])

    # Pad out to max protein length, way more memory than needed but fine here
    max_len = max(len(protein) for protein in proteins_list)
    proteins = np.full((len(proteins_list), max_len), -1, dtype=np.int8)
    lengths = np.empty(len(proteins_list), dtype=np.int32)
    for i, p in enumerate(proteins_list):
        proteins[i, : len(p)] = p
        lengths[i] = len(p)

    pairs = np.asarray(pair_list, dtype=np.int32)

    with ProgressBar(total=len(pairs), unit="pairs") as progress:
        scores = _jitted_score_protein_pairs(
            proteins,
            lengths,
            pairs,
            scoring_matrix,
            config.gap_opening_penalty,
            config.gap_extension_penalty,
            progress,
        )

    # OrthoMCL-style normalization
    self_scores = np.empty(len(proteins_list), dtype=np.float64)
    for i, p in enumerate(proteins_list):
        residues = np.asarray(p, dtype=np.int32)
        self_scores[i] = scoring_matrix[residues, residues].sum()

    return {
        pair: float(scores[i] / np.sqrt(self_scores[pair[0]] * self_scores[pair[1]]))
        for i, pair in enumerate(pair_list)
    }
