"""Cluster perturbation latent vectors and build the Plotly heatmap."""

import argparse
import pathlib

import anndata
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from numba import njit, prange
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from sklearn.metrics.pairwise import cosine_similarity
from statsmodels.stats.multitest import multipletests

from helpers.aggregation import aggregate_anndata


ROOT = pathlib.Path(__file__).resolve().parents[1]
WORK = ROOT / "pipeline_work"


@njit(parallel=True)
def mean_pairwise_distance(x, y):
    total = 0.0
    for i in prange(x.shape[0]):
        for j in range(y.shape[0]):
            squared = 0.0
            for k in range(x.shape[1]):
                difference = x[i, k] - y[j, k]
                squared += difference * difference
            total += np.sqrt(squared)
    return total / (x.shape[0] * y.shape[0])


@njit
def energy_distance(a, b):
    return (
        2 * mean_pairwise_distance(a, b)
        - mean_pairwise_distance(a, a) * (a.shape[0] / (a.shape[0] - 1))
        - mean_pairwise_distance(b, b) * (b.shape[0] / (b.shape[0] - 1))
    )


def permutation_test(a, b, permutations, rng):
    observed = energy_distance(a, b)
    combined = np.vstack([a, b])
    exceedances = 0
    for _ in range(permutations):
        shuffled = rng.permutation(combined.shape[0])
        permuted = energy_distance(
            combined[shuffled[: a.shape[0]]], combined[shuffled[a.shape[0] :]]
        )
        exceedances += permuted >= observed
    return observed, (exceedances + 1) / (permutations + 1)


def compute_edists(embeddings, min_cells, permutations):
    """Test each perturbation against matched Non_target cells within each group."""
    rows = []
    rng = np.random.default_rng(42)
    for group in embeddings.obs["predicted_group"].unique():
        group_data = embeddings[embeddings.obs["predicted_group"] == group]
        controls = group_data[group_data.obs["gene_target"] == "Non_target"].X
        controls = np.asarray(controls, dtype=np.float64)
        if controls.shape[0] > 1000:
            controls = controls[rng.choice(controls.shape[0], 1000, replace=False)]
        if controls.shape[0] < 2:
            continue
        group_rows = []
        targets = group_data.obs["gene_target"].unique()
        for target in targets:
            if target == "Non_target":
                continue
            treated = group_data[group_data.obs["gene_target"] == target].X
            treated = np.asarray(treated, dtype=np.float64)
            if treated.shape[0] < min_cells:
                continue
            distance, pvalue = permutation_test(treated, controls, permutations, rng)
            group_rows.append(
                {"predicted_group": group, "gene_target": target,
                 "edist": distance, "pval": pvalue, "n_cells": treated.shape[0]}
            )
        if group_rows:
            adjusted = multipletests([row["pval"] for row in group_rows], method="fdr_bh")[1]
            for row, pvalue_adjusted in zip(group_rows, adjusted):
                row["pval_adj"] = pvalue_adjusted
            rows.extend(group_rows)
    if not rows:
        raise ValueError("No testable perturbations were found in embeddings.h5ad")
    return pd.DataFrame(rows)


def build_effect_vectors(embeddings, edist, min_cells):
    pseudobulk = aggregate_anndata(
        embeddings, agg_columns=["predicted_group", "gene_target"], sep_char=" - "
    )
    pseudobulk.X = pseudobulk.X / pseudobulk.obs[["n_cells"]].to_numpy()
    controls = pseudobulk[pseudobulk.obs["gene_target"] == "Non_target"]
    controls_by_group = {
        group: np.asarray(controls.X[i]).ravel()
        for i, group in enumerate(controls.obs["predicted_group"])
    }
    effect = pseudobulk[pseudobulk.obs["gene_target"] != "Non_target"].copy()
    effect.X = np.asarray(effect.X, dtype=np.float64)
    for i, group in enumerate(effect.obs["predicted_group"]):
        effect.X[i] -= controls_by_group[group]
    effect = effect[effect.obs["n_cells"] >= min_cells].copy()

    edist = edist.copy()
    edist.index = edist["predicted_group"].astype(str) + " - " + edist["gene_target"].astype(str)
    effect = effect[list(effect.obs_names.intersection(edist.index))].copy()
    effect.obs["edist"] = edist.loc[effect.obs_names, "edist"].to_numpy()
    effect.obs["edist_pval_adj"] = edist.loc[effect.obs_names, "pval_adj"].to_numpy()
    return effect


def write_clustered_html(effect, output_html, significance):
    selected = effect[effect.obs["edist_pval_adj"] < significance]
    if selected.n_obs < 2:
        raise ValueError(f"Only {selected.n_obs} pairs passed FDR < {significance}; clustering needs at least 2")
    labels = selected.obs_names.to_list()
    similarity = cosine_similarity(np.asarray(selected.X, dtype=np.float64))
    order = leaves_list(linkage(pdist(similarity), method="average"))
    labels = np.asarray(labels)[order].tolist()
    similarity = similarity[np.ix_(order, order)]
    figure = go.Figure(go.Heatmap(z=similarity[::-1], x=labels, y=labels[::-1], zmin=-1, zmax=1, colorscale="RdBu_r"))
    figure.update_layout(title="Clustered Perturbation Effect Vectors", width=900, height=900)
    figure.write_html(output_html, include_plotlyjs=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--output-edist", type=pathlib.Path, default=WORK / "edist_results.csv.gz")
    parser.add_argument("--output-effect-vectors", type=pathlib.Path, default=WORK / "pb_effect_vectors.h5ad")
    parser.add_argument("--output-html", type=pathlib.Path, default=WORK / "clustered_pert_effect_vectors.html")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--fdr", type=float, default=0.05)
    args = parser.parse_args()
    for output in (
        args.output_edist,
        args.output_effect_vectors,
        args.output_html,
    ):
        output.parent.mkdir(parents=True, exist_ok=True)

    embeddings = anndata.read_h5ad(args.input)
    required = {"predicted_group", "gene_target"}
    missing = required - set(embeddings.obs.columns)
    if missing:
        raise ValueError(f"embeddings.h5ad is missing .obs columns: {sorted(missing)}")
    edist = compute_edists(embeddings, args.min_cells, args.permutations)
    edist.to_csv(args.output_edist, index=False)
    effect = build_effect_vectors(embeddings, edist, args.min_cells)
    effect.write_h5ad(args.output_effect_vectors)
    write_clustered_html(effect, args.output_html, args.fdr)
    print(args.output_edist.resolve())
    print(args.output_effect_vectors.resolve())
    print(args.output_html.resolve())


if __name__ == "__main__":
    main()
