#!/usr/bin/env python3
# Genome map of ICE community 99 with the prophage community 74 region (Fig 5.7).

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.lines import Line2D

G = "../genomad_genes_for_ch5/"
ICE_LEN, PHAGE_LEN, PHAGE_OFF = 105022, 47313, 66206

def genes(path, key, offset=0):
    out = []
    for i, l in enumerate(open(path)):
        if i == 0:
            continue
        p = l.rstrip("\n").split("\t")
        if len(p) < 20 or key not in p[0]:
            continue
        cls = p[8].split(".")[-1] if p[8] != "NA" else "-"
        out.append(dict(s=int(p[1]) - offset, e=int(p[2]) - offset,
                        cls=cls, desc=("" if p[19] == "NA" else p[19])))
    return out

ice = genes(G + "1307.2545_provirus_genes.tsv", "provirus_114_105135")
phg = genes(G + "1307.3942_provirus_genes.tsv", "provirus_66207_113519", PHAGE_OFF)
hits = []
for l in open("pc74_vs_ice99_hits.tsv"):
    f = l.split()
    hits.append(dict(qs=int(f[0]), qe=int(f[1]), ss=int(f[2]), se=int(f[3])))

PLAS, VIR, NONE, EDGE = "#d98a3d", "#3b6ea5", "#e2e2e2", "#5a5a5a"
col = lambda c: PLAS if c == "PV" else (VIR if c in ("VV", "VP", "VC") else NONE)
rc = lambda x: PHAGE_LEN - x + 1

fig, ax = plt.subplots(figsize=(7.4, 3.5))
BAR, Y_ICE, Y_PHG = 0.115, 0.72, 0.20

for y, ln, gs, flip in ((Y_ICE, ICE_LEN, ice, False), (Y_PHG, PHAGE_LEN, phg, True)):
    ax.add_patch(Rectangle((0, y), ln, BAR, facecolor="white",
                           edgecolor=EDGE, linewidth=0.7, zorder=2))
    for g in gs:
        a, b = (rc(g["e"]), rc(g["s"])) if flip else (g["s"], g["e"])
        ax.add_patch(Rectangle((a, y + 0.004), b - a, BAR - 0.008,
                               facecolor=col(g["cls"]), edgecolor="none", zorder=3))

for h in hits:
    ss, se = sorted((h["ss"], h["se"]))
    qs, qe = sorted((rc(h["qs"]), rc(h["qe"])))
    ax.add_patch(Polygon([(ss, Y_ICE), (se, Y_ICE), (qe, Y_PHG + BAR), (qs, Y_PHG + BAR)],
                         closed=True, facecolor=VIR, alpha=0.13, edgecolor="none", zorder=1))

def bracket(x0, x1, y, txt, style="normal"):
    ax.annotate("", xy=(x0, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle="|-|,widthA=0.3,widthB=0.3",
                                color=EDGE, linewidth=0.7))
    ax.text((x0 + x1) / 2, y - 0.035, txt, ha="center", va="top",
            fontsize=7.2, color=EDGE, style=style)

pv = [g for g in ice if g["cls"] == "PV"]
vr = [g for g in ice if g["cls"] in ("VV", "VP")]
bracket(min(g["s"] for g in pv), max(g["e"] for g in pv), Y_ICE - 0.05,
        "conjugation / plasmid markers", "italic")
bracket(min(g["s"] for g in vr), max(g["e"] for g in vr), Y_ICE - 0.05,
        "phage markers", "italic")

struct = [g for g in ice if "terminase" in g["desc"].lower()
          or "head-tail" in g["desc"].lower() or "tail tube" in g["desc"].lower()
          or "minor tail" in g["desc"].lower() or "gp10" in g["desc"].lower()
          or "gp8" in g["desc"].lower() or "adaptor" in g["desc"].lower()]
s0, s1 = min(g["s"] for g in struct), max(g["e"] for g in struct)
ax.plot([s0, s1], [Y_ICE + BAR + 0.04] * 2, color=EDGE, lw=1.0)
ax.text((s0 + s1) / 2, Y_ICE + BAR + 0.055,
        f"morphogenesis module ({(s1-s0)/1000:.0f} kb)",
        ha="center", fontsize=7.2, color=EDGE)

ax.text(0, Y_ICE + BAR + 0.145, "ICE_community_99 representative  (105,022 bp)",
        fontsize=8.5, weight="bold")
ax.text(0, Y_PHG - 0.115,
        "Prophage community 74 representative, reverse complement  (47,313 bp)",
        fontsize=8.5, weight="bold", va="top")
ax.text(0, Y_PHG - 0.175,
        "40,273 bp (85.1%) aligns to the ICE at 89-95% identity",
        fontsize=7.2, color=EDGE, va="top")

ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=7, color=PLAS,
                          label="plasmid marker"),
                   Line2D([], [], marker="s", ls="", ms=7, color=VIR,
                          label="virus marker"),
                   Line2D([], [], marker="s", ls="", ms=7, color=NONE,
                          label="no marker")],
          loc="lower right", frameon=False, fontsize=7.2,
          handletextpad=0.4, borderpad=0.1, bbox_to_anchor=(1.0, -0.12))

sb = 10000
ax.plot([ICE_LEN - sb, ICE_LEN], [Y_PHG - 0.27] * 2, color=EDGE, lw=1.1)
ax.text(ICE_LEN - sb / 2, Y_PHG - 0.315, "10 kb", ha="center", va="top",
        fontsize=7.2, color=EDGE)

ax.set_xlim(-2500, ICE_LEN + 2500)
ax.set_ylim(-0.22, 1.02)
ax.axis("off")
fig.tight_layout(pad=0.3)
for ext in ("pdf", "png"):
    fig.savefig(f"ICE99_PC74_genome_map.{ext}", dpi=400, bbox_inches="tight",
                facecolor="white")
print("wrote ICE99_PC74_genome_map.pdf / .png")
