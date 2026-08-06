"""Generate per-cell perturbation embeddings from a cleaned raw-count H5AD."""

import argparse
import pathlib

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from scipy.sparse import issparse
from scvi.model import LinearSCVI, SCVI


def require_columns(adata, columns):
    missing = set(columns) - set(adata.obs.columns)
    if missing:
        raise ValueError(f"Input .obs is missing required columns: {sorted(missing)}")


def check_raw_counts(adata):
    """Reject normalized input by checking a small sample of nonzero values."""
    sample = adata[: min(1000, adata.n_obs)].X
    if issparse(sample):
        values = sample.data
    else:
        values = np.asarray(sample).ravel()
        values = values[values != 0]
    if values.size and not np.allclose(values, np.round(values), atol=1e-6):
        raise ValueError("Input .X does not appear to contain raw integer counts")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("embeddings.h5ad"))
    parser.add_argument("--work-dir", type=pathlib.Path, default=pathlib.Path("embedding_models"))
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--latent-dim", type=int, default=100)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument("--scvi-epochs", type=int, default=50)
    parser.add_argument("--linear-scvi-epochs", type=int, default=150)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    adata = anndata.read_h5ad(args.input)
    require_columns(adata, [args.sample_key, "gene_target", "predicted_group"])
    check_raw_counts(adata)

    # Stage 1: standard scVI followed by Leiden clustering.
    SCVI.setup_anndata(adata, batch_key=args.sample_key)
    model = SCVI(adata, n_latent=args.latent_dim, gene_likelihood="nb")
    model.train(
        max_epochs=args.scvi_epochs,
        early_stopping=True,
        early_stopping_warmup_epochs=5,
        batch_size=args.batch_size,
        plan_kwargs={
            "n_epochs_kl_warmup": 5,
            "reduce_lr_on_plateau": True,
            "lr_patience": 4,
            "lr": 0.01,
        },
    )
    model.save(args.work_dir / "scvi_model", overwrite=True)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(
        adata,
        n_neighbors=args.neighbors,
        n_pcs=args.latent_dim,
        use_rep="X_scVI",
    )
    sc.tl.leiden(
        adata,
        resolution=args.leiden_resolution,
        key_added="cluster",
        random_state=args.seed,
    )

    # Stage 2: reproduce the perturbation embedding model using cluster as
    # the categorical covariate. Restrict to single-guide cells when available.
    if "num_guides" in adata.obs:
        linear_ad = adata[adata.obs["num_guides"] == 1].copy()
    else:
        linear_ad = adata.copy()
    require_columns(linear_ad, ["cluster"])
    for internal_column in ["_scvi_batch", "_scvi_labels"]:
        if internal_column in linear_ad.obs:
            del linear_ad.obs[internal_column]
    LinearSCVI.setup_anndata(linear_ad, batch_key="cluster")
    linear_model = LinearSCVI(
        linear_ad,
        n_latent=args.latent_dim,
        gene_likelihood="nb",
    )
    linear_model.train(
        max_epochs=args.linear_scvi_epochs,
        early_stopping=True,
        early_stopping_warmup_epochs=15,
        train_size=0.95,
        batch_size=args.batch_size,
        plan_kwargs={
            "n_epochs_kl_warmup": 15,
            "reduce_lr_on_plateau": True,
            "lr_patience": 10,
            "lr": 0.01,
        },
    )
    linear_model.save(args.work_dir / "linear_scvi_model", overwrite=True)
    latent = linear_model.get_latent_representation()

    embeddings = anndata.AnnData(
        X=latent,
        obs=linear_ad.obs.copy(),
        var=pd.DataFrame(
            index=[f"scVI-{i + 1}" for i in range(latent.shape[1])]
        ),
    )
    embeddings.write_h5ad(args.output)
    print(args.output.resolve())


if __name__ == "__main__":
    main()
