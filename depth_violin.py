#!/usr/bin/env python
# depth_violin.py

import argparse
import scanpy as sc
from scipy import sparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import logging


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def find_genes(adata, sample, region, area, depth_col):
    logging.info(
        f"Finding genes for {sample}–{region}–{area} using depth `{depth_col}`")
    mask = (
        (adata.obs['sample'] == sample) &
        (adata.obs['region'] == region) &
        (adata.obs['area'] == area) &
        (~adata.obs[depth_col].isna()) &
        (adata.obs['H1_annotation'].isin(['EN-ET', 'EN-IT', 'EN-Mig']))
    )
    sub = adata[mask].copy()
    logging.info(f"  Subset: {sub.n_obs} cells × {sub.n_vars} genes")

    X = sub.X.toarray() if sparse.issparse(sub.X) else sub.X
    df = pd.DataFrame({
        'gene': np.repeat(sub.var_names, sub.n_obs),
        'depth': np.tile(sub.obs[depth_col], sub.n_vars)
    })
    ordered = list(df.groupby('gene')['depth'].median().sort_values().index)
    logging.info(f"  Ordered {len(ordered)} genes by median depth")
    return ordered


def make_violin(adata, sample, region, area, genes, depth_col, out_dir):
    logging.info(f"Making violin for {sample}–{region}–{area}")
    mask = (
        (adata.obs['sample'] == sample) &
        (adata.obs['region'] == region) &
        (adata.obs['area'] == area)
    )
    sub = adata[mask].copy()
    if depth_col not in sub.obs.columns or sub.obs[depth_col].isna().all():
        logging.warning(
            f"Depth column `{depth_col}` missing or all-NA for this slice; plotting zeros")
        depths = np.zeros((sub.n_obs,))
    else:
        depths = sub.obs[depth_col].values

    # build data for each gene
    X = sub[:, genes].X
    X = X.toarray() if sparse.issparse(X) else X
    data = [X[:, i] for i in range(X.shape[1])]

    # create violin plot
    fig, ax = plt.subplots(figsize=(len(genes)*0.3, 4))
    sns.violinplot(data=data, inner="quartile", ax=ax)
    ax.set_xticklabels(genes, rotation=90, fontsize=6)
    ax.set_ylabel(depth_col)
    ax.set_title(f"{sample} • {region} • {area}")
    fig.tight_layout()

    fname = f"{out_dir}/{sample}_{region}_{area}_violin.png"
    fig.savefig(fname, dpi=300)
    plt.close(fig)
    logging.info(f"  Saved → {fname}")


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Generate cortical-depth violin plots"
    )
    parser.add_argument('--h5ad', required=True, help="Input AnnData (.h5ad)")
    parser.add_argument('--triples', nargs=3, action='append',
                        metavar=('SAMPLE', 'REGION', 'AREA'),
                        help="One or more SAMPLE REGION AREA triples")
    parser.add_argument('--genes', nargs='+',
                        help="Explicit gene list; skips auto-selection")
    parser.add_argument('--depth-col', default='cortical_depth',
                        help="Column in `adata.obs` to use for depth (fallbacks to zeros)")
    parser.add_argument('--output', default='.', help="Output directory")
    args = parser.parse_args()
    logging.info(f"Arguments: {args}")

    # Load data
    logging.info(f"Loading AnnData from `{args.h5ad}`")
    adata = sc.read(args.h5ad)
    logging.info(f"Available obs columns: {adata.obs.columns.tolist()}")

    # Determine genes
    if args.genes:
        genes = args.genes
        logging.info(f"Using provided gene list ({len(genes)} genes)")
    else:
        s, r, a = args.triples[0]
        try:
            genes = find_genes(adata, s, r, a, args.depth_col)
        except Exception as e:
            logging.error(f"Error while auto-finding genes: {e}")
            sys.exit(1)

    # Plot each triple
    for s, r, a in args.triples:
        make_violin(adata, s, r, a, genes, args.depth_col, args.output)


if __name__ == '__main__':
    main()
