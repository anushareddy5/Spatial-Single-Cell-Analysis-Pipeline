# param_Fig3B_violin.py
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
    log(f"Entering find_genes(sample={sample}, region={region}, area={area})")
    mask = (
        (adata.obs['sample'] == sample) &
        (adata.obs.region == region) &
        (adata.obs.area == area) &
        (~adata.obs['cortical_depth'].isna()) &
        (adata.obs.H1_annotation.isin(['EN-ET', 'EN-IT', 'EN-Mig']))
    )
    adata1 = adata[mask].copy()
    log(f"  Subset size: {adata1.n_obs} cells × {adata1.n_vars} genes")

    X = adata1.X.toarray() if sparse.issparse(adata1.X) else adata1.X
    log("  Computed rounded integer counts")
    df1 = pd.DataFrame({
        'gene': np.repeat(adata1.var_names, adata1.n_obs),
        'cortical_depth': np.tile(adata1.obs['cortical_depth'], adata1.n_vars)
    })
    genes2 = list(
        df1.groupby('gene')['cortical_depth']
           .median()
           .sort_values()
           .index
    )
    log(f"  Final ordered genes: {genes2}")
    return genes2

def make_violin(adata, sample, region, area, genes, out_dir):
    log(f"Entering make_violin(sample={sample}, region={region}, area={area})")
    mask = (
        (adata.obs['sample'] == sample) &
        (adata.obs.region == region) &
        (adata.obs.area == area) &
        (~adata.obs['cortical_depth'].isna()) &
        (adata.obs.H1_annotation.isin(['EN-ET', 'EN-IT', 'EN-Mig']))
    )
    adata1 = adata[mask][:, genes].copy()
    log(f"  Subset for violin: {adata1.n_obs} cells × {len(genes)} genes")

    X = adata1.X.toarray() if sparse.issparse(adata1.X) else adata1.X
    dist_all = []
    for g in genes:
        vals = X[:, adata1.var_names.get_loc(g)]
        dist_all.append(vals)
    fig, ax = plt.subplots(figsize=(len(genes)*0.3, 4))
    sns.violinplot(data=dist_all, inner="quartile", ax=ax)
    ax.set_xticklabels(genes, rotation=90, fontsize=6)
    ax.set_title(f"{sample}-{region}-{area}")
    outfile = f"{out_dir}/{sample}_{region}_{area}_violin.png"
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)
    log(f"  Saved violin to {outfile}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h5ad', required=True)
    p.add_argument('--triples', nargs=3, action='append',
                   metavar=('SAMPLE','REGION','AREA'),
                   help="One or more sample-region-area triples")
    p.add_argument('--genes', nargs='+', help="List of genes (optional)")
    p.add_argument('--output', default='.', help="Output directory")
    args = p.parse_args()

    log("Loading AnnData")
    adata = sc.read(args.h5ad)
    log(f"AnnData loaded: {adata.n_obs}×{adata.n_vars}")

    # If no genes provided, derive them from the first triple
    if not args.genes:
        sample, region, area = args.triples[0]
        genes = find_genes(adata, sample, region, area)
    else:
        genes = args.genes

    for sample, region, area in args.triples:
        make_violin(adata, sample, region, area, genes, args.output)

if __name__ == "__main__":
    main()
