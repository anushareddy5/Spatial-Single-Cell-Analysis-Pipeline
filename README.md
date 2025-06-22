# Spatial Single-Cell Analysis Pipeline

This repository provides a fully parameterized pipeline to reproduce the analysis from the Spatial Single-cell Analysis project using any set of `.h5ad` files, sample–region–area triples, and gene lists.

## Overview & Lineage

The pipeline consists of four core analysis scripts plus a master driver. They must be run in the following order for each `.h5ad` file:

1. **Violin Plot Generation** (`param_Fig3B_violin.py`)
2. **Spatial Gene Expression Plots** (`param_Fig3CDE_spatial_plots_genes.py`)
3. **Expression Heatmaps** (`param_Fig3DE_heatmaps.py`)
4. **Anterior–Posterior Spatial Plots** (`param_Fig3G_spatial_plots_ap.py`)
5. **Master Orchestration** (`run_fig3_pipeline.py`)

Each script reads a common AnnData object, subsets it by sample and region (and area for violins), and then produces publication-quality figures. The master script calls each of the four analysis steps in sequence, ensuring consistent parameter usage and reproducible output organization.

---

## Prerequisites

* Python 3.8+
* Install dependencies:

  ```bash
  pip install scanpy matplotlib seaborn matplotlib-scalebar
  ```
* A directory containing your `.h5ad` files
* A writable output directory

---

## Script Descriptions

### 1. `param_Fig3B_violin.py`

**Purpose**: Generates violin plots of gene expression across cortical depth for selected cell populations.

**Key Features**:

* Subsets cells by `sample`, `region`, and `area` (e.g., `FB123`, `F1`, `A-PFC`).
* Automatically derives an ordered gene list based on median cortical depth if none is provided.
* Saves one violin plot per sample–region–area triple.

**Usage**:

```bash
python param_Fig3B_violin.py \
  --h5ad path/to/data.h5ad \
  --triples SAMPLE REGION AREA \
  [--triples ...] \
  [--genes GENE1 GENE2 ...] \
  [--output output_dir]
```

Arguments:

* `--h5ad`: Path to input `.h5ad` file
* `--triples`: One or more `SAMPLE REGION AREA` triples
* `--genes`: (Optional) List of genes; if omitted, genes are computed for the first triple
* `--output`: (Optional) Output directory (default: current folder)

---

### 2. `param_Fig3CDE_spatial_plots_genes.py`

**Purpose**: Creates spatial gene expression plots for a list of genes across multiple sample–region combinations.

**Key Features**:

* Reads the same `.h5ad` and subsets by `sample` and `region`.
* Scales expression values to a fixed range for comparability.
* Parallelizes plotting across genes with the `--workers` flag.
* Adds a scale bar and color legend.

**Usage**:

```bash
python param_Fig3CDE_spatial_plots_genes.py \
  --h5ad path/to/data.h5ad \
  --samples SAMPLE-REGION [SAMPLE-REGION ...] \
  --genes GENE1 GENE2 ... \
  [--cmap Colormap] \
  [--workers N] \
  [--output output_dir]
```

Arguments:

* `--samples`: List of `SAMPLE-REGION` identifiers (e.g. `FB123-F1`)
* `--genes`: List of genes to plot
* `--cmap`: (Optional) Matplotlib colormap (default: `YlGnBu`)
* `--workers`: (Optional) Number of parallel processes (default: 8)

---

### 3. `param_Fig3DE_heatmaps.py`

**Purpose**: Produces heatmaps of scaled gene expression for each sample–region pair.

**Key Features**:

* Subsets by `sample` and `region`.
* Scales values to a uniform range.
* Arranges genes along the y-axis and cells along the x-axis.

**Usage**:

```bash
python param_Fig3DE_heatmaps.py \
  --h5ad path/to/data.h5ad \
  --samples SAMPLE-REGION [SAMPLE-REGION ...] \
  --genes GENE1 GENE2 ... \
  [--output output_dir]
```

---

### 4. `param_Fig3G_spatial_plots_ap.py`

**Purpose**: Generates anterior–posterior spatial plots (`sc.pl.spatial`) for each gene in each sample–region.

**Key Features**:

* Simple loop over `sample`, `region`, and `genes`.
* Uses Scanpy’s `spatial` plotting function for Visium-style layouts.

**Usage**:

```bash
python param_Fig3G_spatial_plots_ap.py \
  --h5ad path/to/data.h5ad \
  --samples SAMPLE-REGION [SAMPLE-REGION ...] \
  --genes GENE1 GENE2 ... \
  [--output output_dir]
```

---

### 5. `run_fig3_pipeline.py`

**Purpose**: Orchestrates the full pipeline by invoking the four scripts in sequence for each provided `.h5ad`.

**Key Features**:

* Accepts multiple `.h5ad` paths.
* Passes consistent parameters (triples, sample-regions, gene lists) to each step.
* Creates a unified `output` directory structure.

**Usage**:

```bash
python run_fig3_pipeline.py \
  --h5ads file1.h5ad file2.h5ad file3.h5ad \
  --triples FB123 F1 A-PFC \
  --triples O2 Occi A-OC \
  --sr FB123-F1 FB123-O2 \
  --genes CBLN2 SRM RASGRF2 ... \
  --out fig3_results
```

This single command will run:

1. Violin plots for each triple
2. Spatial plots for each `--sr`
3. Heatmaps for each `--sr`
4. AP plots for each `--sr`

All outputs are written under `fig3_results/` with filenames reflecting sample, region, and gene.

---

## Tips & Troubleshooting

* Ensure that the `obs` columns `sample`, `region`, `area`, `cortical_depth`, and `H1_annotation` match exactly between your new `.h5ad` files and the code’s filters.
* If you see memory spikes in spatial plotting, reduce `--workers` in `param_Fig3CDE_spatial_plots_genes.py`.
* For custom gene ordering, always provide `--genes` to skip the automatic median-depth computation.


