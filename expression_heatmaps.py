# expression_heatmaps.py

import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def make_heatmap(adata, sample, region, genes, out_dir):
    sub = adata[(adata.obs['sample'] == sample) & (
        adata.obs['region'] == region)].copy()
    sc.pp.scale(sub, zero_center=True, max_value=6)
    mat = sub[:, genes].X.toarray() if hasattr(sub.X, 'toarray') else sub.X
    df = pd.DataFrame(mat, columns=genes).T

    fig, ax = plt.subplots(figsize=(len(genes)*0.2, 6))
    sns.heatmap(df, cbar=True, ax=ax)
    fig.tight_layout()

    fname = f"{out_dir}/{sample}_{region}_heatmap.png"
    fig.savefig(fname, dpi=300)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h5ad', required=True)
    p.add_argument('--samples', nargs='+', help='SAMPLE-REGION')
    p.add_argument('--genes', nargs='+', required=True)
    p.add_argument('--output', default='.')
    args = p.parse_args()

    adata = sc.read(args.h5ad)
    for sr in args.samples:
        sample, region = sr.split('-')
        make_heatmap(adata, sample, region, args.genes, args.output)


if __name__ == '__main__':
    main()
