#!/usr/bin/env python
# pipeline_runner.py

import argparse
import logging
import subprocess
import sys
import os


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def call(script, args, use_r=False):
    cmd = (['Rscript', script] if use_r else [sys.executable, script]) + args
    logging.info(f"Executing: {' '.join(cmd)}")
    subprocess.check_call(cmd)


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Run full spatial single-cell analysis pipeline"
    )
    parser.add_argument('--h5ads', nargs='+', required=True)
    parser.add_argument('--triples', nargs=3,
                        action='append', metavar=('S', 'R', 'A'))
    parser.add_argument('--sr', nargs='+', required=True,
                        help='List of SAMPLE-REGION strings')
    parser.add_argument('--genes', nargs='+', required=True)
    parser.add_argument('--depth-col', default='cortical_depth',
                        help="Column name in adata.obs to use for depth-violin plots")
    parser.add_argument('--out', default='results')
    parser.add_argument('--spot-size', type=float, default=1.0,
                        help="Spot size for spatial AP plots")
    args = parser.parse_args()

    # Create output dirs
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(f"{args.out}/DEG", exist_ok=True)

    for h in args.h5ads:
        logging.info(f"Processing file: {h}")

        call('depth_violin.py', [
            '--h5ad', h,
            *sum((['--triples', *t] for t in args.triples), []),
            '--genes', *args.genes,
            '--depth-col', args.depth_col,
            '--output', args.out
        ])

        call('spatial_expression_plots.py', [
            '--h5ad', h,
            '--samples', *args.sr,
            '--genes', *args.genes,
            '--spot-size', str(args.spot_size),
            '--output', args.out
        ])

        call('expression_heatmaps.py', [
            '--h5ad', h,
            '--samples', *args.sr,
            '--genes', *args.genes,
            '--output', args.out
        ])

        call('ap_spatial_plots.py', [
            '--h5ad', h,
            '--samples', *args.sr,
            '--genes', *args.genes,
            '--output', args.out
        ])

        call('area_gene_analysis.py', [
            '--h5ad', h,
            '--outdir', f"{args.out}/DEG"
        ])

        # If you have a combined summary CSV, uncomment:
        # call('bubble_plot.R', [
        #     '--input', f"{args.out}/DEG/summary.csv",
        #     '--output', f"{args.out}/DEG/bubble_plot.pdf"
        # ], use_r=True)

    logging.info("All pipeline steps completed successfully.")


if __name__ == '__main__':
    main()
