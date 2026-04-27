from config import PreprocessingConfig, Stage1Config
from preprocessing import preprocess_dataset
from stage_1.kmer_filter import kmer_filter

if __name__ == "__main__":
    preprocessing_config = PreprocessingConfig()
    preprocessed_dataset = preprocess_dataset(preprocessing_config)

    stage_1_config = Stage1Config()
    filtered_kmers = kmer_filter(preprocessed_dataset, stage_1_config)
    print(f"Naive: {len(preprocessed_dataset.proteins) ** 2}")
    print(f"Filtered: {len(filtered_kmers)}")
