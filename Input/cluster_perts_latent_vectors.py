"""
Clustering perturbations using the latent perturbation vectors
"""

from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from numba import njit, prange
from sklearn.metrics.pairwise import cosine_similarity
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

from EffectVectorClustering.Input.aggregation import aggregate_anndata

MIN_CELLS = 20

sns.set_context("paper")

# set fonttype for svg
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


# %% Set up input and output paths

# Path to per-cell scvi embeddings
embeddings_path = "Output/batch1_2_3_4/cluster_perts_latent_vectors/embeddings.h5ad"

out_dir = Path("Output/batch1_2_3_4/cluster_perts_latent_vectors/")
out_dir.mkdir(parents=True, exist_ok=True)

edist_path = out_dir / "edist_results.csv.gz"


# %% Run Edist in this space

# Helpers


@njit(parallel=True)
def mean_pairwise_distance(X, Y):
    nX = X.shape[0]
    nY = Y.shape[0]
    d = 0.0
    for i in prange(nX):
        for j in range(nY):
            dist = 0.0
            for k in range(X.shape[1]):
                diff = X[i, k] - Y[j, k]
                dist += diff * diff
            d += np.sqrt(dist)
    return d / (nX * nY)


@njit
def energy_distance(A, B):
    N = A.shape[0]
    M = B.shape[0]

    term1 = 2 * mean_pairwise_distance(A, B)
    term2 = mean_pairwise_distance(A, A) * (N / (N - 1))
    term3 = mean_pairwise_distance(B, B) * (M / (M - 1))
    return term1 - term2 - term3


def permuted_edist_pvalue(A, B, n_perm=1000, random_state=None):
    rng = np.random.default_rng(random_state)
    X = np.vstack([A, B])
    nA = A.shape[0]
    nB = B.shape[0]
    observed = energy_distance(A, B)
    permuted = np.empty(n_perm)
    for i in range(n_perm):
        idx = rng.permutation(nA + nB)
        A_perm = X[idx[:nA]]
        B_perm = X[idx[nA:]]
        permuted[i] = energy_distance(A_perm, B_perm)
    p_value = (np.sum(permuted >= observed) + 1) / (n_perm + 1)
    return observed, permuted, p_value


def compute_target_edists(ad_embeddings):
    N_NTCS_MAX = 1000
    edist_results_all = {}

    for group in ad_embeddings.obs["predicted_group"].unique():
        if group in edist_results_all:
            continue
        print(group)
        ad_sub_group = ad_embeddings[
            ad_embeddings.obs["predicted_group"] == group
        ].copy()

        valid_ntcs = ad_sub_group.obs.loc[lambda x: x["gene_target"] == "Non_target"]
        if valid_ntcs.shape[0] > N_NTCS_MAX:
            valid_ntcs = valid_ntcs.sample(N_NTCS_MAX, replace=False)
        valid_ntcs = valid_ntcs.index

        perts = set(ad_sub_group.obs["gene_target"].unique()) - set(["Non_target"])

        valid_cells = ad_sub_group.obs.loc[
            lambda x: (
                x.gene_target.isin(perts)
                | ((x.gene_target == "Non_target") & (x.index.isin(valid_ntcs)))
            )
        ]

        ad_sub = ad_sub_group[valid_cells.index]

        ad_sub = ad_sub.copy()

        observed_all = {}
        p_value_all = {}
        n_cells_all = {}
        for gt in tqdm(ad_sub.obs["gene_target"].unique()):
            A = ad_sub[ad_sub.obs["gene_target"] == gt].X

            if A.shape[0] < MIN_CELLS:
                continue

            B = ad_sub[ad_sub.obs["gene_target"] == "Non_target"].X
            observed, permuted, p_value = permuted_edist_pvalue(
                A, B, n_perm=1000, random_state=42
            )
            observed_all[gt] = observed
            # permuted_all[gt] = permuted
            p_value_all[gt] = p_value
            n_cells_all[gt] = A.shape[0]

        edist_results = pd.DataFrame.from_dict(
            {"edist": observed_all, "pval": p_value_all, "n_cells": n_cells_all},
            orient="columns",
        )

        if edist_results.shape[0] == 0:
            continue

        edist_results["pval_adj"] = multipletests(
            edist_results["pval"], method="fdr_bh"
        )[1]
        edist_results["predicted_group"] = group
        edist_results.index.name = "gene_target"
        edist_results = edist_results.reset_index()
        edist_results_all[group] = edist_results

    edist_results_df = pd.concat(edist_results_all.values(), ignore_index=True)

    return edist_results_df


# %% Load embeddings and associated metadata


ad_embeddings = anndata.read_h5ad(embeddings_path)


# %% Pseudobulk the embeddings for each gene_target


ad_embeddings_pb = aggregate_anndata(
    ad_embeddings, agg_columns=["predicted_group", "gene_target"], sep_char=" - "
)

# Divide by num_cells to take average vector
ad_embeddings_pb.X = ad_embeddings_pb.X / ad_embeddings_pb.obs[["n_cells"]].values


# %% Normalize by subtracting out NTC embeddings for each group

# Collect the non-target embeddings

ad_embeddings_pb_ntc = ad_embeddings_pb[
    ad_embeddings_pb.obs["gene_target"] == "Non_target"
]
ntc_embedding_map = {}
for i, group in enumerate(ad_embeddings_pb_ntc.obs["predicted_group"]):
    group_embedding = ad_embeddings_pb_ntc.X[i]
    ntc_embedding_map[group] = group_embedding

# Normalize the perturbation embeddings by taking the difference compared to ntc

ad_embeddings_pb_norm = ad_embeddings_pb[
    ad_embeddings_pb.obs["gene_target"] != "Non_target"
].copy()

for i, group in enumerate(ad_embeddings_pb_norm.obs["predicted_group"]):
    group_embedding = ad_embeddings_pb_norm.X[i]
    ntc_embedding = ntc_embedding_map[group]
    norm_embedding = group_embedding - ntc_embedding
    ad_embeddings_pb_norm.X[i] = norm_embedding

ad_embeddings_pb_norm = ad_embeddings_pb_norm[
    ad_embeddings_pb_norm.obs["n_cells"] >= MIN_CELLS
].copy()

# %% Subset to just gene targets with significant edists

if not edist_path.exists():
    edist_results_df = compute_target_edists(ad_embeddings)
    edist_results_df.to_csv(edist_path, index=True)
else:
    edist_results_df = pd.read_csv(edist_path, index_col=0)

edist_results_df.index = (
    edist_results_df["predicted_group"] + " - " + edist_results_df["gene_target"]
)

ad_embeddings_pb_norm.obs["edist_pval_adj"] = edist_results_df.loc[
    ad_embeddings_pb_norm.obs_names, "pval_adj"
]

ad_embeddings_pb_norm_sub = ad_embeddings_pb_norm[
    ad_embeddings_pb_norm.obs["edist_pval_adj"] < 0.05
].copy()

# %% Plot the subset with the Grin-related targets

gts_to_include = [
    "Grin2a",
    "Grin2b",
    "Grin1",
    "Ap2s1",
    "Ap2m1",
    "Gabrg2",
    "Kcnb1",
    # "Cacna1a",
    # "Mef2a",
    # "Mef2c",
]

groups_to_include = [
    "005 L4-5 IT CTX Glut",
    "007 L2-3 IT CTX Glut",
    "027 NP-CT-L6b-OB Glut",
    "110 CNU-HYa HY MM Glut",
    "151 TH Prkcd Grin2c Glut",
    "155 MB Glut",
    "191 MB P MY GABA",
]

plot_data = ad_embeddings_pb_norm_sub[
    ad_embeddings_pb_norm_sub.obs["gene_target"].isin(gts_to_include)
    & ad_embeddings_pb_norm_sub.obs["predicted_group"].isin(groups_to_include)
].to_df()


# Compute cosine similarity between all pairs of rows
cos_sim_matrix = cosine_similarity(plot_data.values)
cos_sim_df = pd.DataFrame(
    cos_sim_matrix, index=plot_data.index, columns=plot_data.index
)

group_to_coarse_group = {
    "005 L4-5 IT CTX Glut": "Cortex",
    "007 L2-3 IT CTX Glut": "Cortex",
    "027 NP-CT-L6b-OB Glut": "Cortex",
    "110 CNU-HYa HY MM Glut": "Thalamus/Hypothalamus",
    "151 TH Prkcd Grin2c Glut": "Thalamus/Hypothalamus",
    "155 MB Glut": "Midbrain",
    "191 MB P MY GABA": "Midbrain",
}


# Map group codes to colors
set1_colors = plt.get_cmap("Set1").colors
coarse_group_color_map = {
    "Cortex": set1_colors[0],  # dark red
    "Thalamus/Hypothalamus": set1_colors[1],  # blue
    "Midbrain": set1_colors[2],  # green
}


def group_color_map(row_name):
    group_code = row_name.rsplit("-", 1)[0].strip()
    coarse_group = group_to_coarse_group[group_code]
    return coarse_group_color_map[coarse_group]


row_colors = pd.DataFrame(
    {"Tissue": cos_sim_df.index.map(group_color_map)}, index=cos_sim_df.index
)

cm = sns.clustermap(
    cos_sim_df,
    row_cluster=True,
    col_cluster=True,
    vmin=-1,
    vmax=1,
    cmap="vlag",
    yticklabels=True,
    xticklabels=False,
    rasterized=False,
    row_colors=row_colors,
)

# Add legend for row colors
plt.sca(cm.ax_heatmap)
plt.yticks(size=5)
legend_patches = [
    Patch(color=color, label=coarse_group)
    for coarse_group, color in coarse_group_color_map.items()
]

cm.ax_heatmap.legend(
    handles=legend_patches,
    loc="lower left",
    bbox_to_anchor=(1, 1),
    frameon=False,
    title="Tissue",
)
plt.savefig(out_dir / "grin_group.pdf", dpi=200)
plt.close()


# %% Make the Grin2a/2b plot again at the guide level
# %% Create the guide-level embedding vectors

ad_embeddings_pb_guide = aggregate_anndata(
    ad_embeddings, agg_columns=["predicted_group", "guide_call"], sep_char=" - "
)

# Divide by num_cells to take average vector
ad_embeddings_pb_guide.X = (
    ad_embeddings_pb_guide.X / ad_embeddings_pb_guide.obs[["n_cells"]].values
)


# Collect the non-target embeddings

ad_embeddings_pb_guide_ntc = ad_embeddings_pb_guide[
    ad_embeddings_pb_guide.obs["gene_target"] == "Non_target"
]

# This just recalculates the same as before - uses the gene-target groupings to get the ntc pseudobulks
# as a convenience (otherwise we need to compute weighted averages across guides again here - result
# would be identical though)
ntc_embedding_map = {}
for i, group in enumerate(ad_embeddings_pb_ntc.obs["predicted_group"]):
    group_embedding = ad_embeddings_pb_ntc.X[i]
    ntc_embedding_map[group] = group_embedding

# Normalize the other embeddings by taking the difference

ad_embeddings_pb_guide_norm = ad_embeddings_pb_guide[
    ad_embeddings_pb_guide.obs["gene_target"] != "Non_target"
].copy()

for i, group in enumerate(ad_embeddings_pb_guide_norm.obs["predicted_group"]):
    group_embedding = ad_embeddings_pb_guide_norm.X[i]
    ntc_embedding = ntc_embedding_map[group]
    norm_embedding = group_embedding - ntc_embedding
    ad_embeddings_pb_guide_norm.X[i] = norm_embedding

# %% Create the grin plot

gts_to_include = [
    "Grin2a",
    "Grin2b",
    "Grin1",
    "Ap2s1",
    "Ap2m1",
    "Gabrg2",
    "Kcnb1",
    # "Cacna1a",
    # "Mef2a",
    # "Mef2c",
]

groups_to_include = [
    "005 L4-5 IT CTX Glut",
    "007 L2-3 IT CTX Glut",
    "027 NP-CT-L6b-OB Glut",
    "110 CNU-HYa HY MM Glut",
    "151 TH Prkcd Grin2c Glut",
    "155 MB Glut",
    "191 MB P MY GABA",
]

original_plot_rows = ad_embeddings_pb_norm_sub.obs.loc[
    lambda x: x["gene_target"].isin(gts_to_include)
    & x["predicted_group"].isin(groups_to_include),
    ["predicted_group", "gene_target"],
].reset_index(names="group-target")

plot_rows = (
    ad_embeddings_pb_guide_norm.obs.reset_index()[
        ["index", "predicted_group", "gene_target"]
    ]
    .merge(original_plot_rows, how="right")["index"]
    .dropna()
)

plot_data = ad_embeddings_pb_guide_norm[plot_rows].to_df()

# Compute cosine similarity between all pairs of rows
cos_sim_matrix = cosine_similarity(plot_data.values)
cos_sim_df = pd.DataFrame(
    cos_sim_matrix, index=plot_data.index, columns=plot_data.index
)

group_to_coarse_group = {
    "005 L4-5 IT CTX Glut": "Cortex",
    "007 L2-3 IT CTX Glut": "Cortex",
    "027 NP-CT-L6b-OB Glut": "Cortex",
    "110 CNU-HYa HY MM Glut": "Thalamus/Hypothalamus",
    "151 TH Prkcd Grin2c Glut": "Thalamus/Hypothalamus",
    "155 MB Glut": "Midbrain",
    "191 MB P MY GABA": "Midbrain",
}


# Map group codes to colors
set1_colors = plt.get_cmap("Set1").colors
coarse_group_color_map = {
    "Cortex": set1_colors[0],  # dark red
    "Thalamus/Hypothalamus": set1_colors[1],  # blue
    "Midbrain": set1_colors[2],  # green
}


def group_color_map(row_name):
    group_code = row_name.rsplit("-", 1)[0].strip()
    coarse_group = group_to_coarse_group[group_code]
    return coarse_group_color_map[coarse_group]


row_colors = pd.DataFrame(
    {"Tissue": cos_sim_df.index.map(group_color_map)}, index=cos_sim_df.index
)

# Optional - manual re-ordering
REORDER = True
manual_group_target_order = [
    "005 L4-5 IT CTX Glut - Grin1",
    "005 L4-5 IT CTX Glut - Ap2s1",
    "027 NP-CT-L6b-OB Glut - Grin2b",
    "005 L4-5 IT CTX Glut - Grin2b",
    "007 L2-3 IT CTX Glut - Grin2b",
    "007 L2-3 IT CTX Glut - Ap2s1",
    "151 TH Prkcd Grin2c Glut - Grin2b",
    "151 TH Prkcd Grin2c Glut - Ap2s1",
    "151 TH Prkcd Grin2c Glut - Ap2m1",
    "151 TH Prkcd Grin2c Glut - Grin1",
    "191 MB P MY GABA - Ap2s1",
    "155 MB Glut - Ap2s1",
    "110 CNU-HYa HY MM Glut - Ap2s1",
    "110 CNU-HYa HY MM Glut - Ap2m1",
    "191 MB P MY GABA - Ap2m1",
    "155 MB Glut - Ap2m1",
    "191 MB P MY GABA - Grin1",
    "155 MB Glut - Grin2b",
    "155 MB Glut - Grin1",
    "110 CNU-HYa HY MM Glut - Grin2b",
    "155 MB Glut - Gabrg2",
    "191 MB P MY GABA - Gabrg2",
    "027 NP-CT-L6b-OB Glut - Grin2a",
    "005 L4-5 IT CTX Glut - Gabrg2",
    "007 L2-3 IT CTX Glut - Kcnb1",
    "005 L4-5 IT CTX Glut - Kcnb1",
    "005 L4-5 IT CTX Glut - Grin2a",
    "007 L2-3 IT CTX Glut - Grin2a",
    "007 L2-3 IT CTX Glut - Gabrg2",
]

plot_meta = ad_embeddings_pb_guide_norm.obs.loc[plot_rows]


def _map_fun(row):
    key = row["predicted_group"] + " - " + row["gene_target"]
    i = manual_group_target_order.index(key)
    return i


plot_meta["gt_sort"] = plot_meta.apply(_map_fun, axis=1)
order = plot_meta.sort_values(["gt_sort", "guide_call"]).index

if REORDER:
    cos_sim_df = cos_sim_df.loc[order, :].loc[:, order]

cm = sns.clustermap(
    cos_sim_df,
    row_cluster=not REORDER,
    col_cluster=not REORDER,
    vmin=-1,
    vmax=1,
    cmap="vlag",
    yticklabels=True,
    xticklabels=False,
    rasterized=False,
    row_colors=row_colors,
)

# Add legend for row colors
plt.sca(cm.ax_heatmap)
plt.yticks(size=5)
legend_patches = [
    Patch(color=color, label=coarse_group)
    for coarse_group, color in coarse_group_color_map.items()
]

cm.ax_heatmap.legend(
    handles=legend_patches,
    loc="lower left",
    bbox_to_anchor=(1, 1),
    frameon=False,
    title="Tissue",
)
plt.savefig(out_dir / "grin_group_guide-level.pdf", dpi=200)
plt.close()
