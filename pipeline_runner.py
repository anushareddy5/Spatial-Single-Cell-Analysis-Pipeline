# pipeline_runner.py

import argparse
import subprocess
import sys
import logging

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def call(script, args):
    """
    Executes a script (Python or R) with the provided arguments.
    """
    cmd = [sys.executable if script.endswith(
        '.py') else 'Rscript', script] + args
    logging.info(f"Executing command: {' '.join(cmd)}")
    subprocess.check_call(cmd)


def main():
    """
    Main function to parse arguments and orchestrate the pipeline execution.
    """
    # Parse command-line arguments
    p = argparse.ArgumentParser(
        description="Run the spatial single-cell analysis pipeline.")
    p.add_argument('--h5ads', nargs='+', required=True,
                   help="List of input .h5ad files.")
    p.add_argument('--triples', nargs=3, action='append', metavar=('S', 'R', 'A'),
                   help="Triples specifying sample, region, and annotation.")
    p.add_argument('--sr', nargs='+', help="List of SAMPLE-REGION pairs.")
    p.add_argument('--genes', nargs='+', required=True,
                   help="List of genes to analyze.")
    p.add_argument('--out', default='results',
                   help="Output directory for results.")
    args = p.parse_args()

    # Create output directories
    import os
    logging.info(f"Creating output directories at {args.out}")
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(f"{args.out}/DEG", exist_ok=True)

    # Process each input .h5ad file
    for h in args.h5ads:
        logging.info(f"Processing file: {h}")

        # Generate depth violin plots
        logging.info("Generating depth violin plots...")
        call('depth_violin.py', [
            '--h5ad', h,
            *sum((['--triples', *t] for t in args.triples), []),
            '--genes', *args.genes,
            '--output', args.out
        ])

        # Generate spatial expression plots
        logging.info("Generating spatial expression plots...")
        call('spatial_expression_plots.py', [
            '--h5ad', h,
            '--samples', *args.sr,
            '--genes', *args.genes,
            '--output', args.out
        ])

        # Generate expression heatmaps
        logging.info("Generating expression heatmaps...")
        call('expression_heatmaps.py', [
            '--h5ad', h,
            '--samples', *args.sr,
            '--genes', *args.genes,
            '--output', args.out
        ])

        # Generate AP spatial plots
        logging.info("Generating AP spatial plots...")
        call('ap_spatial_plots.py', [
            '--h5ad', h,
            '--samples', *args.sr,
            '--genes', *args.genes,
            '--output', args.out
        ])

        # Perform area gene analysis
        logging.info("Performing area gene analysis...")
        call('area_gene_analysis.py', [
            '--h5ad', h,
            '--outdir', f"{args.out}/DEG"
        ])

        # Generate bubble plot
        logging.info("Generating bubble plot...")
        call('bubble_plot.R', [
            '--input', f"{args.out}/DEG/summary.csv",
            '--output', f"{args.out}/DEG/bubble_plot.pdf"
        ])

    logging.info("Pipeline execution completed successfully.")


if __name__ == '__main__':
    main()
