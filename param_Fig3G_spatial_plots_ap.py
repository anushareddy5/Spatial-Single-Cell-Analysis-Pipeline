# param_Fig3G_spatial_plots_ap.py
import argparse
import scanpy as sc
import matplotlib.pyplot as plt

def make_ap_plot(adata, sample, region, gene, out_dir):
    ad = adata[(adata.obs['sample'] == sample) &
               (adata.obs.region == region)].copy()
    sc.pl.spatial(
        ad,
        color=gene,
        show=False,
        title=f"{sample}-{region}-{gene}"
    )
    outfile = f"{out_dir}/{sample}_{region}_{gene}_AP.png"
    plt.savefig(outfile, dpi=300)
    plt.close()

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
        for gene in args.genes:
            make_ap_plot(adata, sample, region, gene, args.output)

if __name__ == "__main__":
    main()
