import pickle
from pathlib import Path

from config import PreprocessingConfig, Stage1Config
from preprocessing import PreprocessedDataset, preprocess_dataset
from stage_1.kmer_filter import kmer_filter
from stage_1.scoring import score_protein_pairs

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

    print("Scoring protein pairs")
    score_protein_pairs(filtered_protein_pairs, preprocessed_dataset, stage_1_config)
