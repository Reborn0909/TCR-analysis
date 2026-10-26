
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

IN  = sys.argv[1] if len(sys.argv) > 1 else "per_cell_scores.csv.gz"
OUT = sys.argv[2] if len(sys.argv) > 2 else "main_figure_module_scores"

df = pd.read_csv(IN)
GROUPS = ["SS", "Pre-sel", "Post-sel"]
PAL = {"SS": "#8FA9B8", "Pre-sel": "#C08A8A", "Post-sel": "#D0A97E"}
PANELS = [("Exhaustion", "Exhaustion"),
          ("Activation_effector_function", "Activation / Effector function"),
          ("Cytotoxicity", "Cytotoxicity")]

plt.rcParams.update({"font.family": "Helvetica", "font.size": 8,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6})
fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9), constrained_layout=True)

def bracket(ax, x1, x2, y, h, text, fs=7.5):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.7, c="black", clip_on=False)
    ax.text((x1+x2)/2, y+h, text, ha="center", va="bottom", fontsize=fs)

for ax, (col, title) in zip(axes, PANELS):
    data = [df.loc[df["group"] == g, col].dropna().values for g in GROUPS]
    parts = ax.violinplot(data, positions=[1,2,3], widths=0.82,
                          showextrema=False, showmedians=False)
    for pc, g in zip(parts["bodies"], GROUPS):
        pc.set_facecolor(PAL[g]); pc.set_alpha(0.85)
        pc.set_edgecolor("#4d4d4d"); pc.set_linewidth(0.4)
    bp = ax.boxplot(data, positions=[1,2,3], widths=0.14, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=0.8),
                    boxprops=dict(lw=0.5), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    for p, g in zip(bp["boxes"], GROUPS):
        p.set_facecolor(PAL[g]); p.set_alpha(0.95)
    ax.set_xticks([1,2,3]); ax.set_xticklabels(GROUPS, fontsize=8)
    ax.set_ylabel("Score", fontsize=8)
    ax.set_title(title, fontsize=8.5, fontweight="bold", pad=4)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    lo, hi = np.nanmin(df[col]), np.nanmax(df[col]); pad = (hi-lo)*0.30
    steps = [(1,2), (2,3), (1,3)]
    for k, (a, b) in enumerate(steps):
        bracket(ax, a, b, hi + pad*(0.35 + k*0.22), pad*0.06, "****")
    ax.set_ylim(lo - pad*0.4, hi + pad*1.05)

for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}.{ext}", dpi=600 if ext == "png" else None)
print(f"  {OUT}.pdf / .png")
