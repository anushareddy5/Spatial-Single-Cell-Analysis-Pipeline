# depth_violin.py

import argparse
import scanpy as sc
from scipy import sparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def log(msg):
    print(f"[{pd.Timestamp.now()}] {msg}")


def find_genes(adata, sample, region, area):
    log(f"Finding genes for {sample}–{region}–{area}")
    mask = (
        (adata.obs['sample'] == sample) &
        (adata.obs['region'] == region) &
        (adata.obs['area'] == area) &
        (~adata.obs['cortical_depth'].isna()) &
        (adata.obs['H1_annotation'].isin(['EN-ET', 'EN-IT', 'EN-Mig']))
    )
    ad = adata[mask].copy()
    X = ad.X.toarray() if sparse.issparse(ad.X) else ad.X
    df = pd.DataFrame({
        'gene': np.repeat(ad.var_names, ad.n_obs),
        'depth': np.tile(ad.obs['cortical_depth'], ad.n_vars)
    })
    ordered = list(df.groupby('gene')['depth'].median().sort_values().index)
    log(f"  {len(ordered)} genes ordered by median depth")
    return ordered


def make_violin(adata, sample, region, area, genes, out_dir):
    log(f"Making violin for {sample}–{region}–{area}")
    mask = (
        (adata.obs['sample'] == sample) &
        (adata.obs['region'] == region) &
        (adata.obs['area'] == area) &
        (~adata.obs['cortical_depth'].isna()) &
        (adata.obs['H1_annotation'].isin(['EN-ET', 'EN-IT', 'EN-Mig']))
    )
    sub = adata[mask][:, genes].copy()
    X = sub.X.toarray() if sparse.issparse(sub.X) else sub.X
    data = [X[:, sub.var_names.get_loc(g)] for g in genes]

    fig, ax = plt.subplots(figsize=(len(genes)*0.3, 4))
    sns.violinplot(data=data, inner="quartile", ax=ax)
    ax.set_xticklabels(genes, rotation=90, fontsize=6)
    ax.set_title(f"{sample} • {region} • {area}")
    fig.tight_layout()

    fname = f"{out_dir}/{sample}_{region}_{area}_violin.png"
    fig.savefig(fname, dpi=300)
    plt.close(fig)
    log(f"  saved → {fname}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h5ad', required=True)
    p.add_argument('--triples', nargs=3, action='append',
                   metavar=('SAMPLE', 'REGION', 'AREA'))
    p.add_argument('--genes', nargs='+')
    p.add_argument('--output', default='.')
    args = p.parse_args()

    log("Loading AnnData")
    adata = sc.read(args.h5ad)
    if args.genes:
        genes = args.genes
    else:
        s, r, a = args.triples[0]
        genes = find_genes(adata, s, r, a)

    for s, r, a in args.triples:
        make_violin(adata, s, r, a, genes, args.output)


if __name__ == '__main__':
    main()
