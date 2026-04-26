from dataclasses import dataclass


@dataclass
class PreprocessingConfig:
    raw_datasets_path: str = "data/raw_datasets"
    preprocessed_data_output_path: str = "data/preprocessed_data.json"
