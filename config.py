from dataclasses import dataclass


@dataclass
class PreprocessingConfig:
    raw_datasets_path: str = "data/raw_datasets"
    preprocessed_data_output_path: str = "data/preprocessed_data.json"


@dataclass
class Stage1Config:
    kmer_filter_k: int = 5
    kmer_filter_threshold: int = 3

    alignment_scoring_matrix: str = "BLOSUM62"
