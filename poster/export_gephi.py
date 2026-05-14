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

n_genomes = len({p.genome_id for p in dataset.proteins})

CATEGORY_SINGLE_COPY_CORE = 0
CATEGORY_OTHER_CORE = 1
CATEGORY_ACCESSORY = 2
CATEGORY_STRAIN_UNIQUE = 3
CATEGORY_UNCLUSTERED = 4

cluster_to_category: dict[int, int] = {}
for cluster_id, members in clusters.items():
    genome_counts: dict[str, int] = {}
    for protein_idx in members:
        gid = dataset.proteins[protein_idx].genome_id
        genome_counts[gid] = genome_counts.get(gid, 0) + 1
    n_genomes_in_cluster = len(genome_counts)
    if n_genomes_in_cluster == 1:
        cluster_to_category[cluster_id] = CATEGORY_STRAIN_UNIQUE
    elif n_genomes_in_cluster < n_genomes:
        cluster_to_category[cluster_id] = CATEGORY_ACCESSORY
    elif len(members) == n_genomes and all(c == 1 for c in genome_counts.values()):
        cluster_to_category[cluster_id] = CATEGORY_SINGLE_COPY_CORE
    else:
        cluster_to_category[cluster_id] = CATEGORY_OTHER_CORE


with open("data/gephi_nodes.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Id", "Genome", "Length", "Cluster", "Color", "CategoryColor"])
    for i, p in enumerate(dataset.proteins):
        cluster_id = protein_to_cluster.get(i, -1)
        color = cluster_id % N_COLORS if cluster_id >= 0 else N_COLORS
        category_color = cluster_to_category.get(cluster_id, CATEGORY_UNCLUSTERED)
        writer.writerow(
            [i, p.genome_id, len(p.protein_str), cluster_id, color, category_color]
        )

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
