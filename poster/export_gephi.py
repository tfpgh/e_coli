import csv
import pickle

from config import PreprocessingConfig, Stage1Config, Stage2Config
from preprocessing import PreprocessedDataset

N_COLORS = 10
EDGE_THRESHOLD = 0.30

preprocessing_config = PreprocessingConfig()
stage_1_config = Stage1Config()
stage_2_config = Stage2Config()

with open(preprocessing_config.preprocessed_data_output_path) as f:
    dataset = PreprocessedDataset.model_validate_json(f.read())

with open(stage_1_config.scored_pairs_output_path, "rb") as f:
    scored_pairs: dict[tuple[int, int], float] = pickle.load(f)

with open(stage_2_config.clusters_output_path, "rb") as f:
    clusters: dict[int, list[int]] = pickle.load(f)

protein_to_cluster: dict[int, int] = {}
for cluster_id, members in clusters.items():
    for protein_idx in members:
        protein_to_cluster[protein_idx] = cluster_id


with open("data/gephi_nodes.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Id", "Genome", "Length", "Cluster", "Color"])
    for i, p in enumerate(dataset.proteins):
        cluster_id = protein_to_cluster.get(i, -1)
        color = cluster_id % N_COLORS if cluster_id >= 0 else N_COLORS
        writer.writerow([i, p.genome_id, len(p.protein_str), cluster_id, color])

n_edges = 0
with open("data/gephi_edges.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Source", "Target", "Weight", "Type"])
    for (i, j), w in scored_pairs.items():
        if w >= EDGE_THRESHOLD:
            writer.writerow([i, j, f"{w:.4f}", "Undirected"])
            n_edges += 1

print(f"Wrote {len(dataset.proteins):,} nodes to data/gephi_nodes.csv")
print(f"Wrote {n_edges:,} edges to data/gephi_edges.csv")
