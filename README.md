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

Visualizes expression trends along cortical depth by generating violin plots for selected cell populations.

```bash
python depth_violin.py \
  --h5ad path/to/data.h5ad \
  --triples SAMPLE REGION AREA [--triples ...] \
  [--genes GENE1 GENE2 ...] \
  [--output output_dir]
```

* **Inputs**: AnnData `.h5ad`, one or more `SAMPLE REGION AREA` triples
* **Optional**: `--genes` list; if omitted, genes are chosen by median-depth ordering
* **Output**: PNG violin plot per triple

---

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

* **Inputs**: `.h5ad`, `SAMPLE-REGION` identifiers (e.g., `FB123-F1`), gene list
* **Options**: `--cmap` (default `YlGnBu`), `--workers` for parallel plotting
* **Output**: One PNG per (sample–region–gene)

---

### 3. `expression_heatmaps.py`

Produces heatmaps of scaled expression for each sample–region pair.

```bash
python expression_heatmaps.py \
  --h5ad path/to/data.h5ad \
  --samples SAMPLE-REGION [SAMPLE-REGION ...] \
  --genes GENE1 GENE2 ... \
  [--output output_dir]
```

* **Scales**: Zero-centered, max-value clipped
* **Plot**: Genes on y-axis, cells on x-axis

---

### 4. `ap_spatial_plots.py`

Generates anterior–posterior spatial plots using Scanpy’s `spatial` for each gene.

```bash
python ap_spatial_plots.py \
  --h5ad path/to/data.h5ad \
  --samples SAMPLE-REGION [SAMPLE-REGION ...] \
  --genes GENE1 GENE2 ... \
  [--output output_dir]
```

---

### 5. `bubble_plot.R`

Creates bubble plots (e.g., module scores or summary metrics) from a CSV input.

```bash
Rscript bubble_plot.R \
  --input path/to/summary.csv \
  --output path/to/bubble_plot.pdf
```

* **Input CSV** must contain columns: gene/module, sample, region, value1 (size), value2 (color)
* **Dependencies**: `ggplot2`, `reshape2`, `optparse`

---

### 6. `pipeline_runner.py`

**Minimal one‑line example:**

```bash
python pipeline_runner.py --h5ads fileA.h5ad fileB.h5ad --triples SAMPLE REGION AREA --sr SAMPLE-REGION --genes GENE1 GENE2 GENE3 --out output_dir
```

Orchestrates the entire workflow by calling the above scripts in order for each `.h5ad`.

```bash
python pipeline_runner.py \
  --h5ads fileA.h5ad fileB.h5ad ... \
  --triples SAMPLE REGION AREA [--triples ...] \
  --sr SAMPLE-REGION [--sr ...] \
  --genes GENE1 GENE2 ... \
  [--out output_dir]
```

This single command will:

1. Generate depth violins
2. Produce spatial feature plots
3. Draw expression heatmaps
4. Create AP spatial plots
5. (Optionally) run bubble-plot in R if summary CSV is provided

Outputs are organized under `output_dir/` with clear filenames for each step.

---

## Tips & Troubleshooting

* Ensure metadata fields (`sample`, `region`, `area`, `cortical_depth`, `H1_annotation`) exist in your `.h5ad`.
* Adjust `--workers` down if memory is limited during spatial plotting.
* Always specify `--genes` to fix gene order and avoid auto-selection.

*End of README*

## References

* Original pipeline: [https://github.com/ShunzhouJiang/Spatial-Single-cell-Analysis-of-Human-Cortical-Layer-and-Area-Specification/tree/main/Fig3](https://github.com/ShunzhouJiang/Spatial-Single-cell-Analysis-of-Human-Cortical-Layer-and-Area-Specification/tree/main/Fig3)

Add reference: [https://github.com/ShunzhouJiang/Spatial-Single-cell-Analysis-of-Human-Cortical-Layer-and-Area-Specification/tree/main/Fig3](https://github.com/ShunzhouJiang/Spatial-Single-cell-Analysis-of-Human-Cortical-Layer-and-Area-Specification/tree/main/Fig3)
