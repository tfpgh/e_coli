import json
from pathlib import Path

from pydantic import BaseModel

from config import PreprocessingConfig


class PreprocessedGenomeMetadata(BaseModel):
    id: str
    organism: str


class PreprocessedProtein(BaseModel):
    genome_id: str
    protein_str: str


class PreprocessedDataset(BaseModel):
    metadata: list[PreprocessedGenomeMetadata]
    proteins: list[PreprocessedProtein]


def preprocess_dataset(config: PreprocessingConfig) -> PreprocessedDataset:
    metadata: list[PreprocessedGenomeMetadata] = []
    proteins: list[PreprocessedProtein] = []

    for dataset_catalog_path in Path(config.raw_datasets_path).rglob(
        "dataset_catalog.json"
    ):
        data_report_path = dataset_catalog_path.with_name("assembly_data_report.jsonl")
        with open(data_report_path) as f:
            data_report = json.load(f)

            genome_id = data_report["assemblyInfo"]["assemblyName"]
            genome_organism = data_report["assemblyInfo"]["biosample"]["description"][
                "organism"
            ]["organismName"]

            metadata.append(
                PreprocessedGenomeMetadata(
                    id=genome_id,
                    organism=genome_organism,
                )
            )

        with open(dataset_catalog_path) as f:
            dataset_catalog = json.load(f)
            for assembly in dataset_catalog["assemblies"]:
                for file in assembly["files"]:
                    if file["fileType"] == "PROTEIN_FASTA":
                        fasta_path = dataset_catalog_path.parent / file["filePath"]
                        break
                else:
                    continue
                break
            else:
                raise Exception(
                    f"Dataset is missing FASTA file: {dataset_catalog_path}"
                )

        with open(fasta_path) as f:
            for protein_fasta_str in f.read().split(">")[1:]:
                proteins.append(
                    PreprocessedProtein(
                        genome_id=genome_id,
                        protein_str="".join(protein_fasta_str.split("\n")[1:]),
                    )
                )

    return PreprocessedDataset(metadata=metadata, proteins=proteins)
