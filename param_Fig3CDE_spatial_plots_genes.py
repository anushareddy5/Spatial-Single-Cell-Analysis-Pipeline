# param_Fig3CDE_spatial_plots_genes.py
import argparse
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib.colors
from itertools import repeat

def make_plot(adata, cmap, sample, region, gene, out_dir):
    ad = adata[(adata.obs['sample'] == sample) &
               (adata.obs.region == region)].copy()
    sc.pp.scale(ad, zero_center=True, max_value=6)
    fig, ax = plt.subplots()
    sc.pl.embedding(
        ad,
        basis="spatial",
        use_raw=False,
        color=gene,
        show=False,
        s=2,
        color_map=cmap,
        alpha=1,
        colorbar_loc=None,
        ax=ax
    )
    norm = matplotlib.colors.Normalize(
        vmin=ad[:, gene].X.min(),
        vmax=ad[:, gene].X.max()
    )
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, orientation='horizontal', pad=0.02,
                 fraction=0.05, label=gene, shrink=0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    scalebar = ScaleBar(dx=1, units="um", fixed_value=500,
                        location='lower right')
    ax.add_artist(scalebar)
    fig.tight_layout()
    outfile = f"{out_dir}/{sample}_{region}_{gene}.png"
    fig.savefig(outfile, dpi=500)
    plt.close(fig)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h5ad', required=True)
    p.add_argument('--samples', nargs='+',
                   help="List of SAMPLE-REGION strings")
    p.add_argument('--genes', nargs='+', required=True)
    p.add_argument('--cmap', default='YlGnBu')
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--output', default='.')
    args = p.parse_args()

    adata = sc.read(args.h5ad)
    for sr in args.samples:
        sample, region = sr.split('-')
        with Pool(args.workers) as pool:
            pool.starmap(
                lambda s, r, g: make_plot(adata, args.cmap, s, r, g, args.output),
                zip(repeat(sample), repeat(region), args.genes)
            )

if __name__ == "__main__":
    main()
