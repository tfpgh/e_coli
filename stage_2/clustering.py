import numpy as np
import scipy.sparse as sp

from config import Stage2Config


def _build_transition_matrix(
    scored_pairs: dict[tuple[int, int], float], n: int
) -> sp.csc_matrix:
    """
    Convert a dictionary of scored pairs into a sparse matrix:
        transition_matrix[i, j] = weight from i to j
    """
    edges = np.array(list(scored_pairs.keys()), dtype=np.int64)
    weights = np.array(list(scored_pairs.values()), dtype=np.float64)

    source_col = edges[:, 0]
    target_row = edges[:, 1]

    transition_matrix = sp.coo_matrix(
        (weights, (target_row, source_col)),
        shape=(n, n),
    ).tocsc()
    transition_matrix = (transition_matrix + transition_matrix.T).tocsc() + sp.eye(
        n, format="csc"
    )

    return transition_matrix


def _column_normalize(transition_matrix: sp.csc_matrix) -> sp.csc_matrix:
    """Scale each column to sum to 1."""
    col_sums = np.asarray(transition_matrix.sum(axis=0)).ravel()
    return (transition_matrix @ sp.diags(1.0 / col_sums)).tocsc()  # pyright: ignore[reportAttributeAccessIssue]


def _max_column_chaos(transition_matrix: sp.csc_matrix) -> float:
    """
    Per-column chaos = max(column) - sum(column^2).

    As this approaches zero over all columns, clustering is converging.
    """
    col_max = np.asarray(transition_matrix.max(axis=0).todense()).ravel()
    col_sum_of_squares = np.asarray(
        transition_matrix.multiply(transition_matrix).sum(axis=0)
    ).ravel()
    return float(np.max(col_max - col_sum_of_squares))


def _extract_clusters(transition_matrix: sp.csc_matrix) -> dict[int, list[int]]:
    M = transition_matrix.tocsr()
    clusters: dict[int, list[int]] = {}
    for row_idx in range(M.shape[0]):  # pyright: ignore[reportOptionalSubscript]
        members = M.indices[M.indptr[row_idx] : M.indptr[row_idx + 1]]
        if len(members) > 0:
            clusters[int(row_idx)] = sorted(int(m) for m in members)
    return clusters


def cluster_proteins(
    scored_pairs: dict[tuple[int, int], float], n: int, config: Stage2Config
) -> dict[int, list[int]]:
    transition_matrix = _build_transition_matrix(scored_pairs, n)
    transition_matrix = _column_normalize(transition_matrix)

    # Main MCL loop
    for iteration in range(config.mcl_max_iters):
        # Expansion
        transition_matrix_csr = transition_matrix.tocsr()
        expanded = transition_matrix_csr
        for _ in range(config.mcl_expansion - 1):
            expanded = expanded @ transition_matrix_csr

        transition_matrix = sp.csc_matrix(expanded)

        # Inflation
        transition_matrix.data **= config.mcl_inflation

        # Prune, keeps matrix sparse
        transition_matrix.data[transition_matrix.data < config.mcl_prune_threshold] = (
            0.0
        )
        transition_matrix.eliminate_zeros()

        transition_matrix = _column_normalize(transition_matrix)

        # Convergence check
        chaos = _max_column_chaos(transition_matrix)
        print(f"Iter {iteration:3d}: nnz={transition_matrix.nnz:>10} chaos={chaos:.2e}")
        if chaos < config.mcl_convergence_threshold:
            break

    return _extract_clusters(transition_matrix)
