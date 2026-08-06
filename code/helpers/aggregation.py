import anndata
import pandas as pd
import numpy as np
from anndata import AnnData
from anndata.experimental import AnnCollection
from scipy.sparse import csr_matrix, issparse
from tqdm import tqdm


def aggregate_anndata(
    data: AnnData | AnnCollection, agg_columns: list[str], sep_char="_"
) -> AnnData:
    """
    Aggregate single-cell data by summing counts for groups defined by `agg_columns`.

    Assumes that all adatas have the same .var

    Parameters:
    ----------
    data : AnnData | AnnCollection
        Input object with `.obs` containing `agg_columns` and `.X` as the count matrix.

    agg_columns : list[str]
        Columns in `data.obs` used to define groups for aggregation.

    Returns:
    -------
    AnnData
        Aggregated data:
        - `.X`: Dense matrix with summed counts for each group.
        - `.obs`: Metadata for groups, one row per unique combination of `agg_columns`.
        - `.var`: Unchanged from the input.
    """
    # Extract the obs DataFrame
    obs = data.obs

    if not isinstance(obs, pd.DataFrame):  # Occurs if we have an AnnCollectionView
        obs = obs.df

    for col in agg_columns:
        if obs[col].isnull().any():
            raise ValueError(f"Column {col} has null values")

    # Create list of all unique combinations of agg_columns as a list of tuples
    unique_combinations = (
        obs[agg_columns].drop_duplicates().itertuples(index=False, name=None)
    )

    # Create dict that maps agg column key -> integer position
    cols_to_index = {comb: idx for idx, comb in enumerate(unique_combinations)}

    output_X = np.zeros((len(cols_to_index), data.n_vars), dtype="float64")
    cell_counts = np.zeros(len(cols_to_index), dtype="int64")

    batch_size = 10000
    for start in tqdm(range(0, data.n_obs, batch_size)):
        end = min(start + batch_size, data.n_obs)
        batch = data[start:end]

        batch_X, batch_counts = _aggregate_subset(batch, agg_columns, cols_to_index)
        output_X += batch_X
        cell_counts += batch_counts

    # Create an output anndata, use _determine_agg_metadata to create the obs
    agg_obs = _determine_agg_metadata(
        obs, agg_columns, cols_to_index, sep_char=sep_char
    )
    agg_obs["n_cells"] = cell_counts

    if isinstance(data, AnnData):
        var = data.var
    else:
        var = data.adatas[0].var

    return AnnData(X=output_X, obs=agg_obs, var=var)


def _aggregate_subset(data_subset, agg_columns, cols_to_index):
    # Create the sparse aggregator matrix
    row_indices = []
    col_indices = []
    data_values = []

    for row_idx, row in enumerate(
        data_subset.obs[agg_columns].itertuples(index=False, name=None)
    ):
        group_idx = cols_to_index[row]
        row_indices.append(group_idx)
        col_indices.append(row_idx)
        data_values.append(1)

    aggregator_matrix = csr_matrix(
        (data_values, (row_indices, col_indices)),
        shape=(len(cols_to_index), data_subset.shape[0]),
        dtype=data_subset.X.dtype,
    )

    # Perform the aggregation using sparse matrix multiplication
    output_matrix = aggregator_matrix @ data_subset.X

    if issparse(output_matrix):
        output_matrix = output_matrix.toarray()

    # Count cells per group
    cell_counts = np.array(aggregator_matrix.sum(axis=1)).flatten().astype("int64")

    return output_matrix, cell_counts


def _determine_agg_metadata(
    obs: pd.DataFrame, agg_columns: list[str], cols_to_index, sep_char
) -> pd.DataFrame:
    """Creates the output obs dataframe

    For each unique combination of agg columns, determine the
        fields that are the same for all matching rows
    """

    cols_to_keep = (
        obs.groupby(agg_columns, observed=True)
        .nunique()
        .max()
        .loc[lambda x: x <= 1]
        .index
    )

    agg_meta = (
        obs[cols_to_keep.union(agg_columns)].drop_duplicates().set_index(agg_columns)
    )

    cols_to_index_rev = {v: k for k, v in cols_to_index.items()}
    output_order = [cols_to_index_rev[i] for i in range(len(cols_to_index))]

    if len(agg_columns) == 1:
        output_order = [x[0] for x in output_order]
        agg_meta = agg_meta.loc[output_order]
    else:
        agg_meta = agg_meta.loc[output_order].reset_index()
        agg_meta.index = [sep_char.join(x) for x in output_order]

    return agg_meta
