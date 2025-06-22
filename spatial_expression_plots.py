#!/usr/bin/env python
# spatial_expression_plots.py

import argparse
import logging
from multiprocessing import Pool
from itertools import repeat

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib_scalebar.scalebar import ScaleBar


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def plot_one(adata, cmap, sample, region, gene, out_dir):
    logging.info(
        f"Plotting spatial expression of '{gene}' in {sample}-{region}")
    sub = adata[(adata.obs['sample'] == sample) & (
        adata.obs['region'] == region)].copy()
    sc.pp.scale(sub, zero_center=True, max_value=6)
    fig, ax = plt.subplots()
    sc.pl.embedding(
        sub,
        basis='spatial',
        color=gene,
        show=False,
        ax=ax,
        s=2,
        color_map=cmap,
        alpha=1,
        colorbar_loc=None
    )
    norm = matplotlib.colors.Normalize(
        vmin=sub[:, gene].X.min(),
        vmax=sub[:, gene].X.max()
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, orientation='horizontal', pad=0.02,
                 fraction=0.05, label=gene, shrink=0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    ax.add_artist(ScaleBar(dx=1, units='um', fixed_value=500,
                           location='lower right'))
    fig.tight_layout()
    outfile = f"{out_dir}/{sample}_{region}_{gene}.png"
    fig.savefig(outfile, dpi=500)
    plt.close(fig)
    logging.info(f"Saved → {outfile}")


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Spatial expression plots per gene"
    )
    parser.add_argument('--h5ad', required=True)
    parser.add_argument('--samples', nargs='+', required=True,
                        help='List of SAMPLE-REGION strings')
    parser.add_argument('--genes', nargs='+', required=True)
    parser.add_argument('--cmap', default='YlGnBu')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--output', default='.')
    args = parser.parse_args()

    logging.info(f"Loading AnnData: {args.h5ad}")
    adata = sc.read(args.h5ad)
    tasks = []
    for sr in args.samples:
        sample, region = sr.split('-', 1)
        for gene in args.genes:
            tasks.append((adata, args.cmap, sample, region, gene, args.output))

    logging.info(f"Starting multiprocessing with {args.workers} workers")
    with Pool(args.workers) as pool:
        pool.starmap(plot_one, tasks)
    logging.info("Spatial plotting complete.")


if __name__ == '__main__':
    main()
