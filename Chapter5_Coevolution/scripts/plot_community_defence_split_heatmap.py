#!/usr/bin/env python3
# =============================================================================
# Chapter-1-style binomial-exact heatmap, adapted to:
#     prophage community  x  host defence subtype
#
# Chapter 1 was system x system (symmetric) -> upper/lower triangle split.
# This matrix is RECTANGULAR (communities vs subtypes), so the "split" is by
# COLOUR within one grid, matching the chapter-1 palette:
#     ocean blue   = association   (enriched / co-occur,  log2(O/E) > 0)
#     pastel mag.  = dissociation  (depleted / exclusion, log2(O/E) < 0)
#     grey         = not significant (FDR-corrected p >= alpha)
#
#   colour intensity = STRENGTH of interaction = |log2(O/E)|   (effect size)
#       O = community hosts carrying the system
#       E = n_hosts * population frequency of the system
#       NB: intensity gradient is scaled to the SHOWN cells only (relative, not absolute)
#   rows + cols hierarchically clustered -> shared signatures form blocks
#   left strip = #bacterial communities the prophage community spans (binned)
#   row/col filter (MIN_SIG_*) counts SIGNIFICANT INTERACTIONS in either direction
#       (enrichments + depletions together), not associations only
#
# Source: community_defence_enrichment.csv (binomial exact, BH-FDR).
# =============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Patch
import scipy.cluster.hierarchy as sch

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
OUT_DIR   = "figures"
ENR_FILE  = os.path.join("data", "community_defence_enrichment.csv")
SUMM_FILE = os.path.join("data", "prophage_community_summary.csv")
OUT_PNG   = os.path.join(OUT_DIR, "community_defence_split_heatmap_10SIG.png")

MIN_SIG_C = 10        # keep prophage communities with >= this many significant interactions (EITHER direction: enrichments + depletions, q<0.05)
MIN_SIG_S = 10        # keep defence subtypes      with >= this many significant interactions (EITHER direction: enrichments + depletions, q<0.05)
CELL      = 0.20     # inches per heatmap cell
ANNOTATE  = None     # None = auto (annotate observed counts only if the grid is small)
GREY      = "#e8e8e8"

# host-span bins (the chapter hypothesis variable) + a discrete sequential palette
SPAN_BINS   = [0, 1, 2, 4, 9, np.inf]
SPAN_LABELS = ["1", "2", "3-4", "5-9", "10+"]
SPAN_COLORS = ["#f7f7f7", "#fdd49e", "#fc8d59", "#d7301f", "#7f0000"]


def sys_label(s):
    t, _, sub = s.partition(":")
    if not sub:
        return t
    if "unresolved" in sub.lower():
        return f"{t} subtype_unresolved"   # keep the type but flag the unresolved subtype
    return sub


def cluster_order(mat):
    if mat.shape[0] < 3:
        return list(range(mat.shape[0]))
    Z = sch.linkage(mat.values, method="average", metric="euclidean")
    return list(sch.leaves_list(Z))


# ----------------------------------------------------------------------------
# 1. Effect size + significance matrices
# ----------------------------------------------------------------------------
enr  = pd.read_csv(ENR_FILE)
summ = pd.read_csv(SUMM_FILE)

enr["O"] = enr["n_hosts_carrying"]
enr["E"] = enr["n_hosts"] * enr["pop_freq"]
enr["log2OE"] = np.log2(enr["O"].clip(lower=0.5) / enr["E"].clip(lower=0.5))

# 'significant' (q<0.05) is TRUE for BOTH directions, so these counts are total
# significant interactions per community/subtype (enrichments + depletions together),
# NOT associations only. A community with 9 enrichments + 3 depletions counts as 12.
sig = enr[enr["significant"]]
keep_c = sig.groupby("pro_comm").size(); keep_c = keep_c[keep_c >= MIN_SIG_C].index
keep_s = sig.groupby("system").size();   keep_s = keep_s[keep_s >= MIN_SIG_S].index
e = enr[enr["pro_comm"].isin(keep_c) & enr["system"].isin(keep_s)].copy()

val  = e.pivot(index="pro_comm", columns="system", values="log2OE")
sigm = e.pivot(index="pro_comm", columns="system", values="significant").fillna(False).astype(bool)
obs  = e.pivot(index="pro_comm", columns="system", values="O")
clustC = val.where(sigm, 0.0).fillna(0.0)           # signed effect, 0 where ns -> drives clustering

ro = cluster_order(clustC); co = cluster_order(clustC.T)
val, sigm, obs = val.iloc[ro, co], sigm.iloc[ro, co], obs.iloc[ro, co]
rows, cols = list(val.index), list(val.columns)
nr, nc = len(rows), len(cols)
print(f"matrix: {nr} prophage communities x {nc} defence subtypes")
if ANNOTATE is None:
    ANNOTATE = max(nr, nc) <= 45

# ----------------------------------------------------------------------------
# 2. Colour maps (chapter-1 palette) + per-direction normalisation
# ----------------------------------------------------------------------------
assoc_cmap  = mcolors.LinearSegmentedColormap.from_list("ocean_blue",     ["#d6eaf8", "#2e86c1", "#1a5276"])
dissoc_cmap = mcolors.LinearSegmentedColormap.from_list("pastel_magenta", ["#f5d5e0", "#c3447a", "#7b2d50"])

# NB: the colour gradient is normalised to the min/max of the SHOWN cells only.
# Changing MIN_SIG_C/S changes which cells are shown, so the same log2(O/E) value
# can map to a different shade between the 'all' and 'top-N' figures. The per-cell
# maths (O, E, log2OE, significance) is identical; only the colour scaling is relative.
V = val.values; Sg = sigm.values
pos = V[(Sg) & (V > 0)]; neg = -V[(Sg) & (V < 0)]
assoc_norm  = mcolors.Normalize(vmin=float(pos.min()) if pos.size else 0, vmax=float(pos.max()) if pos.size else 1)
dissoc_norm = mcolors.Normalize(vmin=float(neg.min()) if neg.size else 0, vmax=float(neg.max()) if neg.size else 1)
print(f"  significant: {int((Sg & (V>0)).sum())} associations (blue), {int((Sg & (V<0)).sum())} dissociations (magenta)")

# ----------------------------------------------------------------------------
# 3. Draw
# ----------------------------------------------------------------------------
span = dict(zip(summ["pro_comm"], summ["n_bac_communities"]))
def span_color(c):
    v = span.get(c, 1)
    b = np.digitize([v], SPAN_BINS, right=True)[0] - 1
    return SPAN_COLORS[int(np.clip(b, 0, len(SPAN_COLORS) - 1))]

fig_w = nc * CELL + 7
fig_h = nr * CELL + 6
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

for i in range(nr):
    y = nr - 1 - i
    # left strip: host-community span (binned)
    ax.add_patch(Rectangle((-2, y), 1, 1, facecolor=span_color(rows[i]), edgecolor="white", lw=0.3))
    for j in range(nc):
        if Sg[i, j]:
            v = V[i, j]
            color = assoc_cmap(assoc_norm(v)) if v >= 0 else dissoc_cmap(dissoc_norm(-v))
        else:
            color = GREY
        ax.add_patch(Rectangle((j, y), 1, 1, facecolor=color, edgecolor="white", lw=0.3))
        if ANNOTATE and Sg[i, j]:
            br = 0.299 * color[0] + 0.587 * color[1] + 0.114 * color[2]
            ax.text(j + 0.5, y + 0.5, f"{int(obs.values[i, j])}", ha="center", va="center",
                    fontsize=6, color="white" if br < 0.55 else "black")

ax.set_xlim(-2.6, nc); ax.set_ylim(0, nr); ax.set_aspect("equal")

# axis labels: systems on top (rotated), communities on left
fs = max(6.5, min(12, 900 / max(nr, nc)))
ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
ax.set_xticks([j + 0.5 for j in range(nc)])
ax.set_xticklabels([sys_label(c) for c in cols], rotation=90, ha="center", fontsize=fs)
ax.set_yticks([nr - 1 - i + 0.5 for i in range(nr)])
ax.set_yticklabels([r.replace("PRO_com_", "com ") for r in rows], fontsize=fs)
ax.text(-1.5, nr + 0.4, "#", ha="center", va="bottom", fontsize=fs + 3, fontweight="bold")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)

# colourbars (chapter-1 style: association right, dissociation bottom)
sm_a = plt.cm.ScalarMappable(cmap=assoc_cmap, norm=assoc_norm); sm_a.set_array([])
cb_a = fig.colorbar(sm_a, ax=ax, fraction=0.022, pad=0.01, location="right")
cb_a.set_label("Association  log₂(O/E)", fontsize=16); cb_a.ax.tick_params(labelsize=13)
sm_d = plt.cm.ScalarMappable(cmap=dissoc_cmap, norm=dissoc_norm); sm_d.set_array([])
cb_d = fig.colorbar(sm_d, ax=ax, fraction=0.022, pad=0.04, location="bottom", orientation="horizontal")
cb_d.set_label("Dissociation  |log₂(O/E)|", fontsize=16); cb_d.ax.tick_params(labelsize=13)

# legend for the host-span strip (anchored just clear of the association colourbar,
# but close to the axes to avoid excess whitespace)
handles = [Patch(facecolor=SPAN_COLORS[k], edgecolor="grey", label=SPAN_LABELS[k]) for k in range(len(SPAN_LABELS))]
ax.legend(handles=handles, title="# bacterial communities\nthe prophage community spans",
          loc="upper left", bbox_to_anchor=(1.001, 1.0), fontsize=13, title_fontsize=14, frameon=True)

ax.set_title("Prophage community × host defence-subtype  (binomial exact, FDR-corrected p<0.05)\n"
             "blue = association (co-occur)   ·   magenta = dissociation (exclusion)   ·   grey = ns\n"
             "colour intensity = strength |log₂(O/E)|",
             fontsize=18, pad=26)

plt.tight_layout()
fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight", facecolor="white")
print("saved", os.path.abspath(OUT_PNG), f"({fig_w:.0f}x{fig_h:.0f} in, annotate={ANNOTATE})")
