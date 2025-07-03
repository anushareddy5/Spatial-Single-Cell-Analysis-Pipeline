#!/usr/bin/env python
"""
Two PANX1 figures

1. PANX1_H1_all5.pdf
   • single solid curve per H1 (RG, IPC, EN-Mig, EN-IT, EN-ET)
2. PANX1_H1_byLobe_all5.pdf
   • same colours split by cortical lobe
     ─ solid   = Frontal
     ─ dotted  = Parietal
     ─ dashed  = Occipital  (only lobe present in current data)

Curves are 0–1 min-max scaled inside each H1 (or H1×lobe).
"""

import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# ── user paths ────────────────────────────────────────────────
gene = "PANX1"
files = {
    "GW18": "gw18_umb1759.h5ad",
    "GW20": "gw20_umb1031.h5ad",
    "GW21": "gw21_umb1932.h5ad",
}

H1_keep  = ["RG", "IPC", "EN-Mig", "EN-IT", "EN-ET"]
H1_key   = "H1_annotation"
area_key = "area"

# map area codes → lobe
lobe_map = {
    "A-V1": "Occipital", "B-V2": "Occipital", "C-V2": "Occipital",
    "OP": "Occipital", "Occi": "Occipital",
    # add frontal / parietal codes as they appear
}
style_for_lobe = {"Frontal": "solid", "Parietal": "dotted", "Occipital": "dashed"}
# ──────────────────────────────────────────────────────────────

# 1 ─ load & aggregate
rows = []
for gw, path in files.items():
    ad = sc.read(path)
    expr = ad[:, gene].X
    expr = expr.toarray().ravel() if hasattr(expr, "toarray") else np.asarray(expr).ravel()
    md = ad.obs[[H1_key, area_key]].copy()
    md["expr"] = expr
    md["GW"]   = gw
    md = md[md[H1_key].isin(H1_keep)]
    md["lobe"] = md[area_key].map(lobe_map)
    md = md.dropna(subset=["lobe"])
    rows.append(md)
df = pd.concat(rows, ignore_index=True)

mean_lobe = (df.groupby(["GW", H1_key, "lobe"], observed=True)["expr"]
               .mean().reset_index())
mean_all  = (mean_lobe.groupby(["GW", H1_key], observed=True)["expr"]
               .mean().reset_index())

# 2 ─ min-max scale
def mm(s): rng=s.max()-s.min(); return (s-s.min())/rng if rng else 0
mean_all["rel"]  = mean_all.groupby(H1_key)["expr"].transform(mm)
mean_lobe["rel"] = mean_lobe.groupby([H1_key,"lobe"])["expr"].transform(mm)

# colours
palette = dict(zip(H1_keep, sns.color_palette("husl", len(H1_keep))))

# ------------------------------------------------------------------
# Figure 1: all lobes averaged
# ------------------------------------------------------------------
plt.figure(figsize=(4.2,3))
for h1 in H1_keep:
    sub = mean_all[mean_all[H1_key]==h1]
    plt.plot(sub["GW"], sub["rel"], color=palette[h1], lw=2, label=h1)

plt.xlabel("Gestational Week")
plt.ylabel("Relative expression")
plt.title("PANX1 – five H1 groups (avg all lobes)")
plt.xticks(list(files.keys()))
plt.ylim(0, 1.05)                              # head-room
plt.tight_layout(rect=[0, 0, 0.9, 0.96])       # top margin for title
plt.legend(frameon=False)
plt.savefig("PANX1_H1_all5.pdf", dpi=300)
plt.close()
print("Saved PANX1_H1_all5.pdf")

# ------------------------------------------------------------------
# Figure 2: split by lobe
# ------------------------------------------------------------------
plt.figure(figsize=(4.2,3))
for h1 in H1_keep:
    for lobe, grp in mean_lobe[mean_lobe[H1_key]==h1].groupby("lobe"):
        plt.plot(grp["GW"], grp["rel"],
                 color=palette[h1],
                 linestyle=style_for_lobe.get(lobe,"solid"),
                 lw=2)

plt.xlabel("Gestational Week")
plt.ylabel("Relative expression")
plt.title("PANX1 – lobes separated (five H1 groups)")
plt.xticks(list(files.keys()))
plt.ylim(0, 1.05)                              # head-room

# legend
handles_colour = [plt.Line2D([],[],color=palette[h], lw=3) for h in H1_keep]
lobes_present  = mean_lobe["lobe"].unique().tolist()
handles_style  = [plt.Line2D([],[],color="k", lw=3,
                             linestyle=style_for_lobe[l]) for l in lobes_present]
plt.legend(handles_colour+handles_style,
           H1_keep + lobes_present,
           bbox_to_anchor=(1.05,0.5), loc="center left",
           frameon=False, fontsize=8)

plt.tight_layout(rect=[0, 0, 0.9, 0.96])
plt.savefig("PANX1_H1_byLobe_all5.pdf", dpi=300)
plt.close()
print("Saved PANX1_H1_byLobe_all5.pdf")
