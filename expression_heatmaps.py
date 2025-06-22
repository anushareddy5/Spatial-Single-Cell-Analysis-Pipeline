#!/usr/bin/env python
# expression_heatmaps.py

import argparse
import logging

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def make_heatmap(adata, sample, region, genes, out_dir):
    logging.info(f"Generating heatmap for {sample}-{region}")
    sub = adata[
        (adata.obs['sample'] == sample) &
        (adata.obs['region'] == region)
    ].copy()
    sc.pp.scale(sub, zero_center=True, max_value=6)
    mat = sub[:, genes].X
    if hasattr(mat, 'toarray'):
        mat = mat.toarray()
    df = pd.DataFrame(mat, columns=genes).T

    fig, ax = plt.subplots(figsize=(len(genes)*0.2, 6))
    sns.heatmap(df, cbar=True, ax=ax)
    fig.tight_layout()

    out_file = f"{out_dir}/{sample}_{region}_heatmap.png"
    fig.savefig(out_file, dpi=300)
    plt.close(fig)
    logging.info(f"Saved → {out_file}")


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Expression heatmaps per sample-region"
    )
    parser.add_argument('--h5ad', required=True, help="Input .h5ad")
    parser.add_argument('--samples', nargs='+', required=True,
                        help="List of SAMPLE-REGION strings")
    parser.add_argument('--genes', nargs='+', required=True,
                        help="Genes to include")
    parser.add_argument('--output', default='.', help="Output directory")
    args = parser.parse_args()

    logging.info(f"Loading AnnData: {args.h5ad}")
    adata = sc.read(args.h5ad)

    for sr in args.samples:
        sample, region = sr.split('-', 1)
        make_heatmap(adata, sample, region, args.genes, args.output)
    logging.info("Heatmap generation complete.")


if __name__ == '__main__':
    main()
