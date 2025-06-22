#!/usr/bin/env python
# ap_spatial_plots.py

import argparse
import logging

import scanpy as sc
import matplotlib.pyplot as plt


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def make_ap(adata, sample, region, gene, out_dir):
    logging.info(f"Creating AP spatial plot for {gene} in {sample}-{region}")
    sub = adata[(adata.obs['sample'] == sample) & (
        adata.obs['region'] == region)].copy()
    sc.pl.spatial(sub, color=gene, show=False,
                  title=f"{sample}-{region}-{gene}")
    outfile = f"{out_dir}/{sample}_{region}_{gene}_AP.png"
    plt.savefig(outfile, dpi=300)
    plt.close()
    logging.info(f"Saved → {outfile}")


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Anterior–Posterior spatial plots per gene"
    )
    parser.add_argument('--h5ad', required=True)
    parser.add_argument('--samples', nargs='+', required=True,
                        help='List of SAMPLE-REGION strings')
    parser.add_argument('--genes', nargs='+', required=True)
    parser.add_argument('--output', default='.')
    args = parser.parse_args()

    logging.info(f"Loading AnnData: {args.h5ad}")
    adata = sc.read(args.h5ad)
    for sr in args.samples:
        sample, region = sr.split('-', 1)
        for gene in args.genes:
            make_ap(adata, sample, region, gene, args.output)
    logging.info("AP plots complete.")


if __name__ == '__main__':
    main()
