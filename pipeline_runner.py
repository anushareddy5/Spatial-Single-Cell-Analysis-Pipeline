# pipeline_runner.py

import argparse
import subprocess
import sys


def call(script, args):
    cmd = [sys.executable if script.endswith(
        '.py') else 'Rscript', script] + args
    subprocess.check_call(cmd)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--h5ads', nargs='+', required=True)
    p.add_argument('--triples', nargs=3, action='append',
                   metavar=('S', 'R', 'A'))
    p.add_argument('--sr', nargs='+', help='SAMPLE-REGION')
    p.add_argument('--genes', nargs='+', required=True)
    p.add_argument('--out', default='results')
    args = p.parse_args()

    # make dirs
    import os
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(f"{args.out}/DEG", exist_ok=True)

    for h in args.h5ads:
        call('depth_violin.py', [
            '--h5ad', h,
            *sum((['--triples', *t] for t in args.triples), []),
            '--genes', *args.genes,
            '--output', args.out
        ])

        call('spatial_expression_plots.py', [
            '--h5ad', h,
            '--samples', *args.sr,
            '--genes', *args.genes,
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

        call('bubble_plot.R', [
            '--input', f"{args.out}/DEG/summary.csv",
            '--output', f"{args.out}/DEG/bubble_plot.pdf"
        ])


if __name__ == '__main__':
    main()
