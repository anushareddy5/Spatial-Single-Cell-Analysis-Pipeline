#!/usr/bin/env python

import scanpy as sc
import numpy as np
from pathlib import Path

# ───── user config ─────────────────────────────────────────────
h5ads  = [
    "gw18_umb1759.h5ad",
    "gw20_umb1031.h5ad",
    "gw21_umb1932.h5ad",
]
gene      = "PANX1"
vmax_list = [2, 3, 4, 5]        # generate plots 0–2, 0–3, …
spot_size = 25
out_dir   = Path("zscore_plots")
out_dir.mkdir(exist_ok=True)
# ───────────────────────────────────────────────────────────────

for f in h5ads:
    ad = sc.read(f)
    if gene not in ad.var_names:
        print(f"⚠️  {gene} not in {f} – skipping"); continue

    expr = ad[:, gene].X
    expr = expr.toarray().ravel() if hasattr(expr, "toarray") else np.asarray(expr).ravel()
    z = (expr - expr.mean()) / expr.std(ddof=0)
    ad.obs[f"{gene}_z_pos"] = np.clip(z, a_min=0, a_max=None)   # keep positives only

    for vmax in vmax_list:
        sc.pl.spatial(
            ad,
            color=f"{gene}_z_pos",
            vmin=0,
            vmax=vmax,
            spot_size=spot_size,
            cmap="YlGnBu",                  # ← NEW: matches counts palette
            show=False,
            title=f"{gene}  z 0–{vmax}",
            save=None,                      # disable automatic Scanpy name
        )
        # Scanpy puts the figure in plt.gcf(); save it manually
        fname = out_dir / f"show_{Path(f).stem}_{gene}_z0to{vmax}.png"
        import matplotlib.pyplot as plt
        plt.tight_layout()
        plt.savefig(fname, dpi=300)
        plt.close()
        print("Saved", fname)

