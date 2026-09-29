#!/usr/bin/env python3
# Gene map and BLASTn alignment of prophage community 74 vs ICE 139 (Fig S5.1).

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.lines import Line2D

G = "../genomad_genes_for_ch5/"
ICE_LEN, PHAGE_LEN, PHAGE_OFF, ICE_OFF = 212056, 47313, 66206, 1152137

def genes(path, key, offset=0):
    out = []
    for i, l in enumerate(open(path)):
        if i == 0:
            continue
        p = l.rstrip("\n").split("\t")
        if len(p) < 20 or key not in p[0]:
            continue
        out.append(dict(s=int(p[1]) - offset, e=int(p[2]) - offset,
                        cls=(p[8].split(".")[-1] if p[8] != "NA" else "-"),
                        desc=("" if p[19] == "NA" else p[19])))
    return out

ice = genes(G + "1307.1539_provirus_genes.tsv", "provirus_1152138_1364193", ICE_OFF)
phg = genes(G + "1307.3942_provirus_genes.tsv", "provirus_66207_113519", PHAGE_OFF)
hits = [dict(zip(("qs", "qe", "ss", "se"), map(int, l.split()[:4])))
        for l in open("pc74_vs_ice139_hits.tsv")]

PLAS, VIR, NONE, EDGE = "#d98a3d", "#3b6ea5", "#e2e2e2", "#5a5a5a"
col = lambda c: PLAS if c == "PV" else (VIR if c in ("VV", "VP", "VC") else NONE)
rc = lambda x: PHAGE_LEN - x + 1

fig, ax = plt.subplots(figsize=(7.4, 3.4))
BAR, Y_ICE, Y_PHG = 0.115, 0.72, 0.20

for y, ln, gs, flip in ((Y_ICE, ICE_LEN, ice, False), (Y_PHG, PHAGE_LEN, phg, True)):
    ax.add_patch(Rectangle((0, y), ln, BAR, facecolor="white",
                           edgecolor=EDGE, linewidth=0.7, zorder=2))
    for g in gs:
        a, b = (rc(g["e"]), rc(g["s"])) if flip else (g["s"], g["e"])
        ax.add_patch(Rectangle((a, y + 0.004), max(b - a, 250), BAR - 0.008,
                               facecolor=col(g["cls"]), edgecolor="none", zorder=3))

for h in hits:
    ss, se = sorted((h["ss"], h["se"]))
    qs, qe = sorted((rc(h["qs"]), rc(h["qe"])))
    ax.add_patch(Polygon([(ss, Y_ICE), (se, Y_ICE), (qe, Y_PHG + BAR), (qs, Y_PHG + BAR)],
                         closed=True, facecolor=VIR, alpha=0.13, edgecolor="none", zorder=1))

lo = min(min(h["ss"], h["se"]) for h in hits)
hi = max(max(h["ss"], h["se"]) for h in hits)
ax.annotate("", xy=(lo, Y_ICE - 0.055), xytext=(hi, Y_ICE - 0.055),
            arrowprops=dict(arrowstyle="|-|,widthA=0.3,widthB=0.3", color=EDGE, linewidth=0.7))
ax.text((lo + hi) / 2, Y_ICE - 0.10, f"aligned region ({(hi-lo)/1000:.0f} kb)",
        ha="center", va="top", fontsize=7.2, color=EDGE, style="italic")

ax.text(0, Y_ICE + BAR + 0.10, "Largest ICE-like element, genome 1307.1539  (212,056 bp)",
        fontsize=8.5, weight="bold")
ax.text(0, Y_PHG - 0.10,
        "Prophage community 74 representative, reverse complement  (47,313 bp)",
        fontsize=8.5, weight="bold", va="top")
ax.text(0, Y_PHG - 0.16,
        "35,332 bp (74.7%) aligns to the element at 87-95% identity, in 9 blocks",
        fontsize=7.2, color=EDGE, va="top")

ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=7, color=PLAS, label="plasmid marker"),
                   Line2D([], [], marker="s", ls="", ms=7, color=VIR, label="virus marker"),
                   Line2D([], [], marker="s", ls="", ms=7, color=NONE, label="no marker")],
          loc="upper left", frameon=False, fontsize=7.2, ncol=3,
          handletextpad=0.4, columnspacing=1.1, bbox_to_anchor=(0.0, 0.055))

sb = 25000
ax.plot([ICE_LEN - sb, ICE_LEN], [Y_PHG - 0.28] * 2, color=EDGE, lw=1.1)
ax.text(ICE_LEN - sb / 2, Y_PHG - 0.325, "25 kb", ha="center", va="top",
        fontsize=7.2, color=EDGE)

ax.set_xlim(-4000, ICE_LEN + 4000)
ax.set_ylim(-0.22, 1.02)
ax.axis("off")
fig.tight_layout(pad=0.3)
for ext in ("pdf", "png"):
    fig.savefig(f"ICE139_PC74_genome_map.{ext}", dpi=400, bbox_inches="tight", facecolor="white")
print("wrote ICE139_PC74_genome_map.pdf / .png")
