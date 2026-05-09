import pickle
from pathlib import Path

from config import PreprocessingConfig, Stage1Config, Stage2Config
from preprocessing import PreprocessedDataset, preprocess_dataset
from stage_1.kmer_filter import kmer_filter
from stage_1.scoring import score_protein_pairs
from stage_2.clustering import cluster_proteins, summarize_clusters

if __name__ == "__main__":
    preprocessing_config = PreprocessingConfig()
    preprocessing_path = Path(preprocessing_config.preprocessed_data_output_path)
    if preprocessing_path.exists():
        print("Loading preprocessed dataset from JSON")
        with open(preprocessing_path) as f:
            preprocessed_dataset = PreprocessedDataset.model_validate_json(f.read())
    else:
        print("No JSON found, preprocessing dataset")
        preprocessed_dataset = preprocess_dataset(preprocessing_config)
        with open(preprocessing_path, "w") as f:
            f.write(preprocessed_dataset.model_dump_json(indent=4))

    stage_1_config = Stage1Config()
    filtered_proteins_path = Path(stage_1_config.filtered_proteins_output_path)
    if filtered_proteins_path.exists():
        print("Loading protein pairs from pickle")
        with open(filtered_proteins_path, "rb") as f:
            filtered_protein_pairs: list[tuple[int, int]] = pickle.load(f)
    else:
        print("No pickle data found, filtering protein pairs")
        filtered_protein_pairs = kmer_filter(preprocessed_dataset, stage_1_config)
        with open(filtered_proteins_path, "wb") as f:
            pickle.dump(filtered_protein_pairs, f)

    print(f"Naive pair count: {len(preprocessed_dataset.proteins) ** 2:,}")
    print(f"Filtered pair count: {len(filtered_protein_pairs):,}")

    scored_pairs_path = Path(stage_1_config.scored_pairs_output_path)
    if scored_pairs_path.exists():
        print("Loading scored pairs from pickle")
        with open(scored_pairs_path, "rb") as f:
            scored_pairs: dict[tuple[int, int], float] = pickle.load(f)
    else:
        print("No pickle data found, scoring protein pairs")
        scored_pairs = score_protein_pairs(
            filtered_protein_pairs, preprocessed_dataset, stage_1_config
        )
        with open(scored_pairs_path, "wb") as f:
            pickle.dump(scored_pairs, f)

    stage_2_config = Stage2Config()
    clusters_path = Path(stage_2_config.clusters_output_path)
    if clusters_path.exists():
        print("Loading clustered proteins from pickle")
        with open(clusters_path, "rb") as f:
            clustered_proteins: dict[int, list[int]] = pickle.load(f)
    else:
        print("No pickle data found, clustering proteins")
        clustered_proteins = cluster_proteins(
            scored_pairs, len(preprocessed_dataset.proteins), stage_2_config
        )
        with open(clusters_path, "wb") as f:
            pickle.dump(clustered_proteins, f)

    summarize_clusters(clustered_proteins, len(preprocessed_dataset.proteins))
