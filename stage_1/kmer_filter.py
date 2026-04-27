from collections import defaultdict

from tqdm import tqdm

from config import Stage1Config
from preprocessing import PreprocessedDataset


def kmer_filter(
    preprocessed_dataset: PreprocessedDataset, config: Stage1Config
) -> list[tuple[int, int]]:
    k = config.kmer_filter_k
    threshold = config.kmer_filter_threshold
    proteins = preprocessed_dataset.proteins
    n = len(proteins)

    kmer_sets: list[set[str]] = []
    kmer_presence: dict[str, list[int]] = defaultdict(list)
    for protein_idx, protein in enumerate(
        tqdm(proteins, desc="Building kmer filter index")
    ):
        protein_str = protein.protein_str
        protein_kmers = set()
        for i in range(len(protein_str) - k + 1):
            protein_kmers.add(protein_str[i : i + k])

        kmer_sets.append(protein_kmers)

        for kmer in protein_kmers:
            kmer_presence[kmer].append(protein_idx)

    filtered_pairs: list[tuple[int, int]] = []
    for i in tqdm(range(n), desc="Scoring pairs"):
        scores: dict[int, int] = defaultdict(int)
        for kmer in kmer_sets[i]:
            shared_kmer_proteins = kmer_presence[kmer]
            for j in shared_kmer_proteins:
                if j > i:
                    scores[j] += 1

        for j, score in scores.items():
            if score >= threshold:
                filtered_pairs.append((i, j))

    return filtered_pairs
