from dataclasses import dataclass


@dataclass
class PreprocessingConfig:
    raw_datasets_path: str = "data/raw_datasets"
    preprocessed_data_output_path: str = "data/preprocessed_data.json"


@dataclass
class Stage1Config:
    kmer_filter_k: int = 5
    kmer_filter_threshold: int = 3
    filtered_proteins_output_path: str = "data/filtered_protein_pairs.pkl"
    scored_pairs_output_path: str = "data/scored_protein_pairs.pkl"

    alignment_scoring_matrix: str = "BLOSUM62"
    gap_opening_penalty: float = 11.0
    gap_extension_penalty: float = 1.0


@dataclass
class Stage2Config:
    mcl_inflation: float = 2.0
    mcl_expansion: int = 2
    mcl_prune_threshold: float = 1e-6
    mcl_edge_threshold: float = 0.20
    mcl_max_iters: int = 100
    mcl_convergence_threshold: float = 1e-6
    clusters_output_path: str = "data/clustered_proteins.pkl"
