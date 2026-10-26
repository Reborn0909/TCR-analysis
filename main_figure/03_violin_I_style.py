
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN  = sys.argv[1] if len(sys.argv) > 1 else "per_cell_scores.csv.gz"
OUT = sys.argv[2] if len(sys.argv) > 2 else "main_figure_module_scores_Istyle"

df = pd.read_csv(IN)
GROUPS = ["SS", "Pre-sel", "Post-sel"]
PAL = {"SS": "#8FA9B8", "Pre-sel": "#C08A8A", "Post-sel": "#D0A97E"}
PANELS = [("Exhaustion", "Exhaustion Signature"),
          ("Activation_effector_function", "Activation / Effector Signature"),
          ("Cytotoxicity", "Cytotoxicity Signature")]
YLO, YHI = -1.0, 2.0          

plt.rcParams.update({
    "font.family": "Helvetica", "font.size": 9,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
})
FIG_W, FIG_H = 3.5, 3.5/(456/708)   
fig, axes = plt.subplots(3, 1, figsize=(FIG_W, FIG_H), sharey=True,
                         constrained_layout=True)

def bracket(ax, x1, x2, y, h, text, fs=9):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.8, c="black", clip_on=False)
    ax.text((x1+x2)/2, y+h, text, ha="center", va="bottom", fontsize=fs)


STEPS = [(1, 2, 1.14), (2, 3, 1.34), (1, 3, 1.56)]   

for ax, (col, title) in zip(axes, PANELS):
    data = [df.loc[df["group"] == g, col].dropna().values for g in GROUPS]
    parts = ax.violinplot(data, positions=[1, 2, 3], widths=0.78,
                          showextrema=False, showmedians=False)
    for pc, g in zip(parts["bodies"], GROUPS):
        pc.set_facecolor(PAL[g]); pc.set_alpha(0.85); pc.set_zorder(3)
        pc.set_edgecolor("#4d4d4d"); pc.set_linewidth(0.45)
    bp = ax.boxplot(data, positions=[1, 2, 3], widths=0.13, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=0.9),
                    boxprops=dict(lw=0.55), whiskerprops=dict(lw=0.55),
                    capprops=dict(lw=0.55))
    for p_, g in zip(bp["boxes"], GROUPS):
        p_.set_facecolor('white'); p_.set_alpha(1.0); p_.set_zorder(4)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(GROUPS, fontsize=9.5)
    ax.set_ylabel("Score", fontsize=9.5)
    ax.set_title(title, fontsize=9.5, fontweight="bold", pad=7)
    ax.grid(True, axis="both", color="#E8E8E8", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.set_xlim(0.45, 3.55)
    ax.set_ylim(YLO, YHI)                       
    ax.set_yticks([-1, 0, 1, 2])                      
    ax.set_yticklabels(["-1", "0", "1", "2"], fontsize=9.5)
    a = 0.5                                     
    for x1, x2, y in STEPS:
        bracket(ax, x1, x2, y, 0.035, "****", fs=9)

for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}.{ext}", dpi=600 if ext == "png" else None)
print(f" → {OUT}.pdf / .png  (figsize {FIG_W:.2f}×{FIG_H:.2f} in =  {FIG_W/FIG_H:.4f}, ylim {YLO}~{YHI})")
