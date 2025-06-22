# spatial_expression_plots.py

import argparse
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from multiprocessing import Pool
from itertools import repeat
import matplotlib.colors


def plot_one(adata, cmap, sample, region, gene, out_dir):
    sub = adata[(adata.obs['sample'] == sample) & (
        adata.obs['region'] == region)].copy()
    sc.pp.scale(sub, zero_center=True, max_value=6)
    fig, ax = plt.subplots()
    sc.pl.embedding(
        sub, basis='spatial', color=gene, show=False, ax=ax,
        s=2, color_map=cmap, alpha=1, colorbar_loc=None
    )
    norm = matplotlib.colors.Normalize(
        vmin=sub[:, gene].X.min(), vmax=sub[:, gene].X.max()
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, orientation='horizontal', pad=0.02,
                 fraction=0.05, label=gene, shrink=0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    ax.add_artist(ScaleBar(dx=1, units='um',
                  fixed_value=500, location='lower right'))
    fig.tight_layout()

    fname = f"{out_dir}/{sample}_{region}_{gene}.png"
    fig.savefig(fname, dpi=500)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h5ad', required=True)
    p.add_argument('--samples', nargs='+', help='SAMPLE-REGION')
    p.add_argument('--genes', nargs='+', required=True)
    p.add_argument('--cmap', default='YlGnBu')
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--output', default='.')
    args = p.parse_args()

    adata = sc.read(args.h5ad)
    tasks = []
    for sr in args.samples:
        sample, region = sr.split('-')
        for g in args.genes:
            tasks.append((adata, args.cmap, sample, region, g, args.output))

    with Pool(args.workers) as pool:
        pool.starmap(plot_one, tasks)


if __name__ == '__main__':
    main()
