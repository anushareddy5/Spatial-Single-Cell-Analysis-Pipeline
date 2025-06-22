#!/usr/bin/env python
# area_gene_analysis.py

from anndata import read_h5ad
import scanpy as sc
import pandas as pd
import numpy as np
import logging
import argparse
import os
os.environ['PYTHONHASHSEED'] = '0'


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def change_name(name):
    return name.split('-', 1)[-1] if '-' in name else name


def expr_tot(adata):
    areas = np.unique(adata.obs['area'])
    logging.info(f"Computing mean expression for areas: {areas}")
    mat = [adata[adata.obs['area'] == a].X.mean(axis=0) for a in areas]
    return pd.DataFrame(np.array(mat).T, index=adata.var.index, columns=areas)


def prop_tot(adata):
    areas = np.unique(adata.obs['area'])
    logging.info(f"Computing detection proportion for areas: {areas}")
    mat = [np.mean(adata[adata.obs['area'] == a].X != 0, axis=0)
           for a in areas]
    return pd.DataFrame(np.array(mat).T, index=adata.var.index, columns=areas)


def run_de(adata, groupby, group, ref, top_n, prefix, expr_df, prop_df):
    logging.info(f"DE: {group} vs {ref} ({groupby})")
    ad = adata.copy()
    sc.tl.rank_genes_groups(ad, groupby, method='t-test',
                            groups=[group], reference=ref)
    names = pd.DataFrame(ad.uns['rank_genes_groups']['names'][:top_n])
    expr_df.loc[names[group], :].to_csv(f"{prefix}_expr_{group}.csv")
    prop_df.loc[names[group], :].to_csv(f"{prefix}_prop_{group}.csv")
    logging.info(f"Saved DE files with prefix '{prefix}'")


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Area-wise expression & DE analysis"
    )
    parser.add_argument('--h5ad', required=True, help="Input .h5ad")
    parser.add_argument('--outdir', default='result/DEG', help="Output folder")
    args = parser.parse_args()

    logging.info(f"Loading AnnData: {args.h5ad}")
    adata = read_h5ad(args.h5ad)
    adata.obs['area'] = adata.obs['area'].map(change_name)
    adata = adata[adata.obs['area'].isin(
        ['PFC', 'V2', 'Par', 'V1', 'M1', 'Temp'])]
    adata = adata[adata.obs['gw'].isin(['gw20', 'gw22'])]
    adata_sub = adata[adata.obs['H1_annotation'].isin(
        ['EN-IT', 'EN-ET'])].copy()

    os.makedirs(args.outdir, exist_ok=True)

    zs_tot = expr_tot(adata_sub)[['PFC', 'M1', 'Par', 'Temp', 'V2', 'V1']]
    pr_tot = prop_tot(adata_sub)[['PFC', 'M1', 'Par', 'Temp', 'V2', 'V1']]

    # A vs P
    adata_sub.obs['direction'] = np.where(
        adata_sub.obs['area'].isin(['PFC', 'M1']), 'A',
        np.where(adata_sub.obs['area'].isin(['Par', 'V2']), 'P', '-1')
    )
    run_de(adata_sub[adata_sub.obs['direction'] != '-1'],
           'direction', 'A', 'P', 50, f"{args.outdir}/AP", zs_tot, pr_tot)

    # Temp vs N
    adata_sub.obs['direction'] = np.where(
        adata_sub.obs['area'] == 'Temp', 'Temp', 'N'
    )
    run_de(adata_sub[adata_sub.obs['direction'] != '-1'],
           'direction', 'Temp', 'N', 20, f"{args.outdir}/Temp", zs_tot, pr_tot)

    # by H1 classes
    for cls in ["EN-Mig", "RG", "IPC", "IN"]:
        sub = adata[adata.obs['H1_annotation'] == cls]
        zs = expr_tot(sub)[['PFC', 'M1', 'Par', 'Temp', 'V2', 'V1']]
        pr = prop_tot(sub)[['PFC', 'M1', 'Par', 'Temp', 'V2', 'V1']]
        zs.to_csv(f"{args.outdir}/expr_{cls}.csv")
        pr.to_csv(f"{args.outdir}/prop_{cls}.csv")
        logging.info(f"Saved class-wise CSVs for {cls}")

    logging.info("Area-gene analysis complete.")


if __name__ == '__main__':
    main()
