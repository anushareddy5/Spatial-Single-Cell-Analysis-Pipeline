# param_Fig3DE_heatmaps.py
import argparse
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

def make_heatmap(adata, sample, region, genes, out_dir):
    ad = adata[(adata.obs['sample'] == sample) &
               (adata.obs.region == region)].copy()
    sc.pp.scale(ad, zero_center=True, max_value=6)
    data = ad[:, genes].X.toarray() if hasattr(ad.X, "toarray") else ad.X
    df = pd.DataFrame(data, columns=genes)
    fig, ax = plt.subplots(figsize=(len(genes)*0.2, 6))
    sns.heatmap(df.T, cbar=True, ax=ax)
    outfile = f"{out_dir}/{sample}_{region}_heatmap.png"
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h5ad', required=True)
    p.add_argument('--samples', nargs='+',
                   help="List of SAMPLE-REGION strings")
    p.add_argument('--genes', nargs='+', required=True)
    p.add_argument('--output', default='.')
    args = p.parse_args()

    adata = sc.read(args.h5ad)
    for sr in args.samples:
        sample, region = sr.split('-')
        make_heatmap(adata, sample, region, args.genes, args.output)

if __name__ == "__main__":
    main()
