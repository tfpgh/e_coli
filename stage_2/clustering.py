import numpy as np
import scipy.sparse as sp

from config import Stage2Config


def _build_transition_matrix(
    scored_pairs: dict[tuple[int, int], float], n: int, edge_threshold: float
) -> sp.csc_matrix:
    """
    Convert a dictionary of scored pairs into a sparse matrix:
        transition_matrix[i, j] = weight from i to j
    """
    edges = np.array(list(scored_pairs.keys()), dtype=np.int64)
    weights = np.array(list(scored_pairs.values()), dtype=np.float64)

    mask = weights >= edge_threshold
    edges, weights = edges[mask], weights[mask]
    print(f"Edges after threshold: {len(weights):,}")

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
    """Scale each column to sum to 1. Zero-sum columns are left as zero."""
    col_sums = np.asarray(transition_matrix.sum(axis=0)).ravel()
    inv = np.zeros_like(col_sums)
    np.divide(1.0, col_sums, out=inv, where=col_sums > 0)
    return (transition_matrix @ sp.diags(inv)).tocsc()  # pyright: ignore[reportAttributeAccessIssue]


def _extract_clusters(transition_matrix: sp.csc_matrix) -> dict[int, list[int]]:
    """Each node is assigned to its dominant attractor (column-wise argmax)."""
    M = transition_matrix.tocsc()
    n_cols = M.shape[1]  # pyright: ignore[reportOptionalSubscript]
    clusters: dict[int, list[int]] = {}
    for col_idx in range(n_cols):
        start, end = M.indptr[col_idx], M.indptr[col_idx + 1]
        if start == end:
            continue
        best_local = M.data[start:end].argmax()
        attractor = int(M.indices[start + best_local])
        clusters.setdefault(attractor, []).append(int(col_idx))
    return clusters


def cluster_proteins(
    scored_pairs: dict[tuple[int, int], float], n: int, config: Stage2Config
) -> dict[int, list[int]]:
    transition_matrix = _build_transition_matrix(
        scored_pairs, n, config.mcl_edge_threshold
    )
    transition_matrix = _column_normalize(transition_matrix)

    prev_matrix: sp.csc_matrix | None = None
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

        if prev_matrix is None:
            max_change = float("inf")
        else:
            diff = transition_matrix - prev_matrix
            max_change = float(np.abs(diff.data).max()) if diff.nnz > 0 else 0.0
        print(
            f"Iter {iteration:3d}: nnz={transition_matrix.nnz:>10} Δmax={max_change:.2e}"
        )
        if max_change < config.mcl_convergence_threshold:
            break

        prev_matrix = transition_matrix.copy()

    return _extract_clusters(transition_matrix)
