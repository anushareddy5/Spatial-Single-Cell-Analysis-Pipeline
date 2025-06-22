# ap_spatial_plots.py

import argparse
import scanpy as sc
import matplotlib.pyplot as plt


def make_ap(adata, sample, region, gene, out_dir):
    sub = adata[(adata.obs['sample'] == sample) & (
        adata.obs['region'] == region)].copy()
    sc.pl.spatial(sub, color=gene, show=False,
                  title=f"{sample}-{region}-{gene}")
    fname = f"{out_dir}/{sample}_{region}_{gene}_AP.png"
    plt.savefig(fname, dpi=300)
    plt.close()


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
        for g in args.genes:
            make_ap(adata, sample, region, g, args.output)


if __name__ == '__main__':
    main()
