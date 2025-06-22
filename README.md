# Spatial Single-Cell Analysis Pipeline

This repository offers a generic, fully parameterized pipeline for spatial single-cell expression analysis. You can apply it to any set of `.h5ad` files, sample–region–area triples, and gene lists without altering core logic.

## Overview & Workflow

Six scripts form the core of this pipeline. They must be executed in the sequence below for each `.h5ad` file:

1. **Depth Violin Plots** (`depth_violin.py`)
2. **Spatial Expression Plots** (`spatial_expression_plots.py`)
3. **Expression Heatmaps** (`expression_heatmaps.py`)
4. **Anterior–Posterior Spatial Plots** (`ap_spatial_plots.py`)
5. **Bubble Plot Analysis** (`bubble_plot.R`)
6. **Pipeline Runner** (`pipeline_runner.py`)

Each step reads the same AnnData object, subsets by the appropriate metadata fields, and produces high-quality visual summaries.

---

## Prerequisites

### Python dependencies

* Python 3.8+
* Install Python libraries:

  ```bash
  pip install scanpy matplotlib seaborn matplotlib-scalebar
  ```

### R dependencies

* R 4.0+

* Install R packages:

  ```r
  install.packages(c("ggplot2", "reshape2", "optparse"))
  ```

* A directory of input `.h5ad` files

* A writable output directory

---

## Script Details

### 1. `depth_violin.py`

**Optional**: `--depth-col` specifies which `adata.obs` column to use for depth (default: `cortical_depth`).

Visualizes expression trends along cortical depth by generating violin plots for selected cell populations.

```bash
python depth_violin.py \
  --h5ad path/to/data.h5ad \
  --triples SAMPLE REGION AREA [--triples ...] \
  [--genes GENE1 GENE2 ...] \
  [--depth-col <depth_column_name>]  # replace with the actual depth column in your data \
  [--output output_dir]
```

* **Inputs**: AnnData `.h5ad`, one or more `SAMPLE REGION AREA` triples
* **Optional**:

  * `--genes`: list of genes; if omitted, genes are auto-selected by median-depth ordering
  * `--depth-col`: obs column name for depth (default `cortical_depth`)
  * `--output`: output directory (default `.`)
* **Output**: PNG violin plot per triple

### 2. `spatial_expression_plots.py`

Creates spatial feature plots for a list of genes across sample–region combinations.

```bash
python spatial_expression_plots.py \
  --h5ad path/to/data.h5ad \
  --samples SAMPLE-REGION [SAMPLE-REGION ...] \
  --genes GENE1 GENE2 ... \
  [--cmap Colormap] \
  [--workers N] \
  [--output output_dir]
```

* **Inputs**: `.h5ad`, `SAMPLE-REGION` identifiers, gene list
* **Options**:

  * `--cmap`: Matplotlib colormap (default `YlGnBu`)
  * `--workers`: number of parallel workers (default 8)
  * `--output`: output directory (default `.`)
* **Output**: PNG per sample–region–gene

### 3. `expression_heatmaps.py`

Produces heatmaps of scaled expression for each sample–region pair.

```bash
python expression_heatmaps.py \
  --h5ad path/to/data.h5ad \
  --samples SAMPLE-REGION [SAMPLE-REGION ...] \
  --genes GENE1 GENE2 ... \
  [--output output_dir]
```

* **Scales**: zero-centered, max-value clipped
* **Plot**: genes on y-axis, cells on x-axis

### 4. `ap_spatial_plots.py`

**Optional**: `--spot-size` sets the point size for spatial plots (default: `1.0`).

Generates anterior–posterior spatial plots using Scanpy’s `spatial` for each gene.

```bash
python ap_spatial_plots.py \
  --h5ad path/to/data.h5ad \
  --samples SAMPLE-REGION [SAMPLE-REGION ...] \
  --genes GENE1 GENE2 ... \
  [--spot-size size] \
  [--output output_dir]
```

### 5. `bubble_plot.R`

Creates bubble plots (e.g., module scores or summary metrics) from a CSV input.

```bash
Rscript bubble_plot.R \
  --input path/to/summary.csv \
  --output path/to/bubble_plot.pdf
```

* **Input CSV** must contain: gene/module, sample, region, value1 (size), value2 (color)
* **Dependencies**: `ggplot2`, `reshape2`, `optparse`

### 6. `pipeline_runner.py`

**Optional**: `--depth-col`, `--spot-size` flags default to `cortical_depth` and `1.0`.

**Minimal one-line example:**

```bash
python pipeline_runner.py \
  --h5ads fileA.h5ad fileB.h5ad \
  --triples SAMPLE REGION AREA \
  --sr SAMPLE-REGION \
  --genes GENE1 GENE2 GENE3 \
  --depth-col <depth_column_name>  # replace with your actual depth-column \
  --spot-size 1.0 \
  --out output_dir
```

Orchestrates the full workflow:

```bash
python pipeline_runner.py \
  --h5ads fileA.h5ad fileB.h5ad fileC.h5ad \
  --triples SAMPLE REGION AREA [--triples ...] \
  --sr SAMPLE-REGION [--sr ...] \
  --genes GENE1 GENE2 ... \
  --depth-col <depth_column_name>  # replace with your actual depth-column \
  --spot-size 1.0 \
  [--out output_dir]
```

This runs:

1. Depth violins
2. Spatial feature plots
3. Expression heatmaps
4. AP spatial plots
5. (Optionally) bubble-plot in R

Outputs saved under `output_dir/`.

---

## Checking Gene Presence and Identifying Sample–Region–Area Triples

To ensure a gene is present in your `.h5ad` and to discover which sample–region–area combinations actually express it, run this Python snippet:

```python
import scanpy as sc
import numpy as np

# 1. Load your AnnData
adata = sc.read("path/to/data.h5ad")

# 2. Specify the gene to check
gene = "PANX1"

# 3. Verify the gene is in var_names
if gene not in adata.var_names:
    print(f"{gene} not found in this dataset.")
else:
    print(f"{gene} found. Identifying expressing cells...")
    # 4. Subset to cells with expression > 0
    expr = adata[:, gene].X
    expr = expr.toarray().flatten() if hasattr(expr, "toarray") else expr.flatten()
    cells_idx = np.where(expr > 0)[0]
    adata_pos = adata[cells_idx].copy()
    # 5. List unique sample-region-area triples
    combos = adata_pos.obs[['sample','region','area']].drop_duplicates()
    print("Sample–Region–Area combinations with", gene, "> 0:")
    print(combos.to_string(index=False))
```

## Tips & Troubleshooting

* Ensure metadata fields (`sample`, `region`, `area`, `cortical_depth`, `H1_annotation`) exist.
* Adjust `--workers` down if memory is limited.
* Always specify `--genes` to fix gene order.
* Check available obs columns with `print(adata.obs.columns)` if depth-col warnings occur.
* You can also list them directly using:

  ```python
  import scanpy as sc
  adata = sc.read("file.h5ad")
  print(adata.obs.columns.tolist())
  ```

## References

* [Original pipeline on GitHub](https://github.com/ShunzhouJiang/Spatial-Single-cell-Analysis-of-Human-Cortical-Layer-and-Area-Specification/tree/main/Fig3)
