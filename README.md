# Effect-vector clustering pipeline

This pipeline converts a cleaned raw-count H5AD into an interactive
low-resolution perturbation-effect module browser.

## Required input

`input.h5ad` must contain raw integer counts in `.X` and the following `.obs`
columns:

- `sample_id`: sample/batch identifier
- `gene_target`: perturbation target; controls must be named `Non_target`
- `predicted_group`: cell-type/group annotation
- `num_guides`: optional; if present, only `num_guides == 1` cells are used

## Run the pipeline step by step

```bash
cd /data/EffectVectorClustering
python -m pip install -r requirements.txt

# Step 1: Train SCVI models and generate per-cell latent embeddings.
python code/generate_embeddings.py \
  --input /path/to/input.h5ad \
  --output pipeline_work/embeddings.h5ad

# Step 2: Calculate perturbation effect vectors and create the clustered heatmap.
python code/cluster_perts_latent_vectors.py \
  --input pipeline_work/embeddings.h5ad \
  --output-edist pipeline_work/edist_results.csv.gz \
  --output-effect-vectors pipeline_work/pb_effect_vectors.h5ad \
  --output-html pipeline_work/clustered_pert_effect_vectors.html

# Step 3: Define low-resolution modules and create the module browser.
python code/generate_low_resolution_browser.py \
  --input pipeline_work/clustered_pert_effect_vectors.html \
  --output-tables analysis_outputs/tables/low_resolution_modules.xlsx \
  --output-html analysis_outputs/html/low_resolution_module_browser.html

# Step 3: Define low-resolution modules and create the module browser.
python code/run_go_enrichment.py \
  --input analysis_outputs/html/low_resolution_module_browser.html \
  --output-tables analysis_outputs/tables/final_low_resolution_modules.xlsx \
  --output-html analysis_outputs/html/final_low_resolution_module_browser.html
```

## Code order

### 1. `generate_embeddings.py`

- **Input:** cleaned `input.h5ad` with raw counts in `.X`
- **Function:** trains standard SCVI using `sample_id`, generates Leiden
  clusters, and trains LinearSCVI using the generated `cluster` annotation
- **Output:** `embeddings.h5ad` containing 100-dimensional per-cell latent
  vectors; trained models under `embedding_models/`

### 2. `cluster_perts_latent_vectors.py`

- **Input:** `embeddings.h5ad`
- **Function:** performs energy-distance permutation tests, averages embeddings
  by `predicted_group × gene_target`, subtracts the matched `Non_target` mean,
  filters significant pairs, calculates cosine similarity, and performs
  average-linkage hierarchical clustering
- **Outputs:**
  - `pipeline_work/edist_results.csv.gz`
  - `pipeline_work/pb_effect_vectors.h5ad`
  - `pipeline_work/clustered_pert_effect_vectors.html`

### 3. `generate_low_resolution_browser.py`

- **Input:** clustered cosine-similarity matrix embedded in
  `clustered_pert_effect_vectors.html`
- **Function:** divides the hierarchical order into contiguous low-resolution
  modules and embeds their summary and member information in the browser
- **Outputs:**
  - `analysis_outputs/tables/low_resolution_modules.xlsx`, with `summary` and
    `members` sheets
  - `analysis_outputs/html/low_resolution_module_browser.html`

### 4. `run_go_enrichment.py`

- **Input:** low-resolution module-browser HTML, which contains the embedded
  module summary and membership data from step 3
- **Function:** queries mouse g:Profiler for `GO:BP`, `GO:MF`, and `GO:CC`,
  using all heatmap perturbation genes as the custom background; adds each
  module's top GO term to the browser
- **Outputs:**
  - `analysis_outputs/tables/final_low_resolution_modules.xlsx`, with `summary`
    and `go_enrichment` sheets
  - `analysis_outputs/html/final_low_resolution_module_browser.html`

## Helper files

The following internal modules reside under `code/helpers/`.

### `helpers/aggregation.py`

- **Input:** per-cell AnnData and grouping columns
- **Function:** calculates pseudobulk sums and cell counts
- **Output:** aggregated AnnData used by step 2

### `helpers/clustering_helpers.py`

- **Input:** clustered Plotly heatmap
- **Function:** reads its similarity matrix, segments it into modules, and
  calculates module statistics
- **Output:** module summaries and membership records used by step 3

### `helpers/module_browser.py`

- **Input:** similarity matrix, module summaries, and module memberships
- **Function:** builds the self-contained interactive HTML interface
- **Output:** low-resolution module-browser HTML used as the final browser

## Etc.

- Energy-distance testing uses 1,000 permutations and may take time on a large dataset.
- GO enrichment requires internet access to the g:Profiler API.
