from config import PreprocessingConfig
from preprocessing import preprocess_dataset

if __name__ == "__main__":
    preprocessing_config = PreprocessingConfig()
    preprocessed_dataset = preprocess_dataset(preprocessing_config)
