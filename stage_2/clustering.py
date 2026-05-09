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


def summarize_clusters(clusters: dict[int, list[int]], n_proteins: int) -> None:
    sizes = np.array([len(members) for members in clusters.values()])
    total_membership = int(sizes.sum())

    print(f"\nClusters: {len(clusters):,}")
    print(f"Total membership: {total_membership:,} (proteins: {n_proteins:,})")
    if total_membership > n_proteins:
        print(f"  → overlap: {total_membership / n_proteins:.2f}x (soft MCL convergence)")
    elif total_membership < n_proteins:
        print(f"  → unassigned: {n_proteins - total_membership:,}")
    print(f"Largest cluster: {int(sizes.max()):,}")
    print(f"Median cluster size: {int(np.median(sizes))}")
    print(f"Singletons: {int((sizes == 1).sum()):,}")
    print("Size distribution:")
    buckets: list[tuple[str, np.ndarray]] = [
        ("= 1", sizes == 1),
        ("2-5", (sizes >= 2) & (sizes <= 5)),
        ("6-20", (sizes >= 6) & (sizes <= 20)),
        ("21-50", (sizes >= 21) & (sizes <= 50)),
        ("51+", sizes >= 51),
    ]
    for label, mask in buckets:
        print(f"  {label:>6}: {int(mask.sum()):,}")


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

        chaos = _max_column_chaos(transition_matrix)
        if prev_matrix is None:
            max_change = float("inf")
        else:
            diff = transition_matrix - prev_matrix
            max_change = float(np.abs(diff.data).max()) if diff.nnz > 0 else 0.0
        print(
            f"Iter {iteration:3d}: nnz={transition_matrix.nnz:>10} "
            f"chaos={chaos:.2e} Δmax={max_change:.2e}"
        )
        if max_change < config.mcl_convergence_threshold:
            break

        prev_matrix = transition_matrix.copy()

    return _extract_clusters(transition_matrix)
