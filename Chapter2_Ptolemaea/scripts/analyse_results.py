#!/usr/bin/env python3
# ======================================================================
# Thesis Chapter 2, Figure 2.1 and Results
# All Chapter 2 statistics (results_stats.txt, e.g. 32,509 annotations, 50.6% multi-source), system counts table, landscape and concordance panels.
# Original location: github.com/ecampbell50/Chp4_Ptolemaea scripts/analyse_results.py
# ======================================================================
"""
analyse_results.py — panel-wide analysis + figures for the Ptolemaea demonstration.

Reads the per-species *_annotations.csv and *_summary.tsv files in the repo root
and produces:
  * results_stats.txt   — every number cited in the Demonstration section
  * figures/fig2_landscape.{pdf,png}
  * figures/fig3_concordance.{pdf,png}

Usage:  python3 scripts/analyse_results.py
"""
import os
import textwrap
from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns  # registers 'rocket_r' etc. as matplotlib colormaps
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# species file-tag -> display name, in a sensible (ESKAPE + E. coli) order
SPECIES = OrderedDict([
    ("Ecoli",        "E. coli"),
    ("Efaecium",     "E. faecium"),
    ("Saureus",      "S. aureus"),
    ("Kpneumoniae",  "K. pneumoniae"),
    ("Abaumanii",    "A. baumannii"),
    ("Paeruginosa",  "P. aeruginosa"),
    ("Enterobacter", "Enterobacter spp."),
])

NULL = "No_hit"
CURATED_STATUSES = {"MAPPING", "CONFLICT"}

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
def clean_type(s):
    """Harmonise a known curation typo: CRISP-Cas -> CRISPR-Cas."""
    return s.replace("CRISP-Cas", "CRISPR-Cas") if isinstance(s, str) else s

ann = {}
summ = {}
for tag in SPECIES:
    a = pd.read_csv(os.path.join(ROOT, "data", f"{tag}_annotations.csv"), dtype=str).fillna(NULL)
    a["final_type"] = a["final_type"].map(clean_type)
    a["species"] = SPECIES[tag]
    ann[tag] = a
    s = pd.read_csv(os.path.join(ROOT, "data", f"{tag}_summary.tsv"), sep="\t")
    s["species"] = SPECIES[tag]
    summ[tag] = s

ALL = pd.concat(ann.values(), ignore_index=True)
SUM = pd.concat(summ.values(), ignore_index=True)

# source membership (PADLOC / DefenseFinder / BLAST) per protein
ALL["has_padloc"] = ALL["padloc_original"] != NULL
ALL["has_df"]     = ALL["deffind_original"] != NULL
ALL["has_blast"]  = (ALL["fwd_blast"] != NULL) | (ALL["rev_blast"] != NULL)
ALL["n_sources"]  = ALL[["has_padloc", "has_df", "has_blast"]].sum(axis=1)

NAMED_EXCLUDE = {"type_unresolved", "Other"}  # catch-all / unresolved buckets

# ---------------------------------------------------------------------------
# Headline numbers
# ---------------------------------------------------------------------------
n_genomes = SUM.shape[0]
total_hits = len(ALL)
distinct_types_all = ALL["final_type"].nunique()
distinct_types_named = ALL.loc[~ALL["final_type"].isin(NAMED_EXCLUDE), "final_type"].nunique()
distinct_subtypes = ALL.loc[~ALL["final_subtype"].isin(
    {"subtype_unresolved", "type_unresolved"}), "final_subtype"].nunique()

status_counts = ALL["status"].value_counts()
status_pct = (status_counts / total_hits * 100).round(1)

multi = (ALL["n_sources"] >= 2).sum()
multi_pct = round(multi / total_hits * 100, 1)

# per-genome systems
pg = SUM["n_defence_genes"]
pg_types = SUM["n_unique_types"]

# recall: how many proteins each single source recovers vs the union
n_padloc = int(ALL["has_padloc"].sum())
n_df     = int(ALL["has_df"].sum())
n_blast  = int(ALL["has_blast"].sum())
best_single = max(n_padloc, n_df, n_blast)
best_name = {"padloc": n_padloc, "DefenseFinder": n_df, "BLAST": n_blast}
best_label = max(best_name, key=best_name.get)
union = total_hits  # every retained protein has >=1 source
extra_over_best = union - best_single
extra_pct = round(extra_over_best / best_single * 100, 1)

# top named systems (exclude catch-all/unresolved)
top_named = (ALL.loc[~ALL["final_type"].isin(NAMED_EXCLUDE), "final_type"]
             .value_counts().head(15))

# curated concentration
curated = ALL[ALL["status"].isin(CURATED_STATUSES)]
curated_top = curated["final_type"].value_counts().head(10)

# Venn regions (PADLOC=P, DefenseFinder=D, BLAST=B)
def region(p, d, b):
    return ((ALL["has_padloc"] == p) & (ALL["has_df"] == d) & (ALL["has_blast"] == b)).sum()
regions = {
    "P_only":  region(True, False, False),
    "D_only":  region(False, True, False),
    "B_only":  region(False, False, True),
    "PD":      region(True, True, False),
    "PB":      region(True, False, True),
    "DB":      region(False, True, True),
    "PDB":     region(True, True, True),
}

# ---------------------------------------------------------------------------
# Write stats report
# ---------------------------------------------------------------------------
lines = []
def w(s=""):
    lines.append(s)

w("=" * 70)
w("PTOLEMAEA DEMONSTRATION — PANEL-WIDE STATISTICS")
w("=" * 70)
w(f"Genomes analysed:                 {n_genomes}")
w(f"Total defence-system annotations: {total_hits:,}")
w(f"Distinct system types (final_type, incl. Other/unresolved): {distinct_types_all}")
w(f"Distinct *named* system types (excl. Other & type_unresolved): {distinct_types_named}")
w(f"Distinct subtypes (excl. unresolved): {distinct_subtypes}")
w(f"Multi-source support (>=2 of PADLOC/DF/BLAST): {multi:,} ({multi_pct}%)")
w("")
w("Status composition (each protein has exactly one status):")
for st in ["AGREE", "RESOLVED", "SINGLE", "BLAST", "MAPPING", "CONFLICT"]:
    c = int(status_counts.get(st, 0))
    w(f"   {st:<9} {c:>7,}  ({status_pct.get(st, 0):>4}%)")
w(f"   {'curated (MAPPING+CONFLICT)':<9} "
  f"{int(status_counts.get('MAPPING',0)+status_counts.get('CONFLICT',0)):>7,}  "
  f"({round((status_counts.get('MAPPING',0)+status_counts.get('CONFLICT',0))/total_hits*100,1)}%)")
w("")
w("Defence systems per genome (n_defence_genes):")
w(f"   min={pg.min()}  max={pg.max()}  median={pg.median():.0f}  mean={pg.mean():.1f}")
w("Distinct system TYPES per genome (n_unique_types):")
w(f"   min={pg_types.min()}  max={pg_types.max()}  median={pg_types.median():.0f}  mean={pg_types.mean():.1f}")
w("")
w("Per species — median (min-max) defence genes / genome:")
for tag, disp in SPECIES.items():
    sub = SUM[SUM["species"] == disp]["n_defence_genes"]
    w(f"   {disp:<18} median={sub.median():.0f}  ({sub.min()}-{sub.max()})  mean={sub.mean():.1f}")
w("")
w("RECALL ARGUMENT (proteins recovered by each source):")
w(f"   PADLOC:        {n_padloc:,}")
w(f"   DefenseFinder: {n_df:,}")
w(f"   BLAST (bidir): {n_blast:,}")
w(f"   Union (all 3): {union:,}")
w(f"   Best single source = {best_label} ({best_single:,})")
w(f"   -> Union recovers {extra_over_best:,} more annotations than the best "
  f"single source alone (+{extra_pct}%).")
w("")
w("Source-membership regions (for Venn, Fig 3A):")
for k, v in regions.items():
    w(f"   {k:<8} {v:>7,}  ({round(v/total_hits*100,1)}%)")
w(f"   (check sum = {sum(regions.values()):,} == total {total_hits:,})")
w("")
w("Top 15 NAMED system types (panel-wide count):")
for name, c in top_named.items():
    w(f"   {name:<16} {c:>6,}")
w("")
w("Where curated (MAPPING/CONFLICT) entries concentrate (top 10 types):")
for name, c in curated_top.items():
    w(f"   {name:<16} {c:>6,}")
w("=" * 70)

report = "\n".join(lines)
print(report)
with open(os.path.join(ROOT, "results_stats.txt"), "w") as fh:
    fh.write(report + "\n")

# ---------------------------------------------------------------------------
# Supplementary Table S2 — per-system abundance across species
# ---------------------------------------------------------------------------
disp_cols = list(SPECIES.values())
s2 = (ALL.groupby(["final_type", "final_subtype", "final_outcome", "species"])
      .size().unstack("species", fill_value=0))
for c in disp_cols:
    if c not in s2.columns:
        s2[c] = 0
s2 = s2[disp_cols]
s2["Total"] = s2.sum(axis=1)
# number of genomes in which each system appears
ngen = (ALL.groupby(["final_type", "final_subtype", "final_outcome"])["genome_id"]
        .nunique().rename("n_genomes"))
s2 = s2.join(ngen)
s2 = s2.sort_values("Total", ascending=False).reset_index()
s2 = s2.rename(columns={"final_type": "type", "final_subtype": "subtype",
                        "final_outcome": "outcome"})
s2.to_csv(os.path.join(ROOT, "SupplementaryTableS2_system_counts.tsv"),
          sep="\t", index=False)
print(f"\nSupplementary Table S2: {len(s2)} system rows "
      f"-> SupplementaryTableS2_system_counts.tsv")

# ===========================================================================
# FIGURES
# ===========================================================================
plt.rcParams.update({
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,   # editable text in PDF
    "ps.fonttype": 42,
})
disp_order = list(SPECIES.values())
palette = dict(zip(disp_order, plt.cm.tab10(np.linspace(0, 1, 10))[:len(disp_order)]))

# ---------------------------------------------------------------------------
# FIGURE 2 — landscape, compared across three sources
#   columns: full Ptolemaea consensus | PADLOC alone | DefenseFinder alone
#   row 1 (A-C): defence genes per genome by species (box + jitter)
#   row 2 (D-F): heatmap of top consensus SUBTYPES (mean copies / genome)
# All panels share axes/colour scales so the three sources are comparable.
# Per-tool panels use the harmonised consensus subtype names, restricted to the
# proteins that tool detected (padloc_original / deffind_original != No_hit).
# Cells with a mean of 0 (subtype completely absent in that species) are drawn
# white, so a true miss is visually distinct from a low-but-present count.
# ---------------------------------------------------------------------------
import numpy.ma as ma
from matplotlib.colors import LinearSegmentedColormap

sp_of = SUM.set_index("genome_id")["species"]
genome_order = sp_of.index
n_per_sp = {sp: int((sp_of == sp).sum()) for sp in disp_order}

# rows: most abundant named subtypes (exclude unresolved buckets + 'Other' catch-all)
SUBTYPE_EXCLUDE = {"subtype_unresolved", "type_unresolved"}
TOP_N = 22
top_subtypes = (ALL.loc[~ALL["final_subtype"].isin(SUBTYPE_EXCLUDE)
                        & (ALL["final_type"] != "Other"), "final_subtype"]
                .value_counts().head(TOP_N).index.tolist())

def per_genome_counts(mask):
    return ALL.loc[mask].groupby("genome_id").size().reindex(genome_order, fill_value=0)

def heat_matrix(mask):
    sub = ALL.loc[mask]
    m = np.zeros((len(top_subtypes), len(disp_order)))
    for j, sp in enumerate(disp_order):
        s = sub[sub["species"] == sp]
        vc = s[s["final_subtype"].isin(top_subtypes)]["final_subtype"].value_counts()
        for i, t in enumerate(top_subtypes):
            m[i, j] = vc.get(t, 0) / n_per_sp[sp]
    return m

# viridis colour scale; 0 cells (absent) are masked and drawn white
CMAP = matplotlib.colormaps["viridis"].copy()
CMAP.set_bad("white")

# (title, second line, protein mask) — explicit about what each column shows
SOURCES = [
    ("Ptolemaea consensus", "(PADLOC + DefenseFinder + BLAST)", np.ones(len(ALL), dtype=bool)),
    ("BLAST-only",          "(genes detected using BLASTp)",    ALL["has_blast"].values),
    ("PADLOC only",         "(genes detected by PADLOC)",       ALL["has_padloc"].values),
    ("DefenseFinder only",  "(genes detected by DefenseFinder)", ALL["has_df"].values),
]
pg = {name: per_genome_counts(mask) for name, _, mask in SOURCES}
hm = {name: heat_matrix(mask) for name, _, mask in SOURCES}

print("\nPer-genome defence genes by source (median [min-max]):")
for name, _, _m in SOURCES:
    c = pg[name]
    print(f"   {name:<22} median={c.median():.0f}  ({c.min()}-{c.max()})  total={int(c.sum()):,}")

ymax = max(s.max() for s in pg.values())
sqrt_vmax = max(np.sqrt(m).max() for m in hm.values())

fig2, axes = plt.subplots(2, 4, figsize=(20, 10.8),
                          gridspec_kw={"height_ratios": [1.0, 2.25]},
                          constrained_layout=True)
letters = [["A", "B", "C", "D"], ["E", "F", "G", "H"]]
im = None
for col, (name, subtitle, _) in enumerate(SOURCES):
    # ---- row 1: per-genome distribution ----
    axA = axes[0, col]
    data_by_sp = [pg[name][sp_of == sp].values for sp in disp_order]
    bp = axA.boxplot(data_by_sp, vert=True, widths=0.62, patch_artist=True,
                     showfliers=False, medianprops=dict(color="black", lw=1.3))
    for patch, sp in zip(bp["boxes"], disp_order):
        patch.set_facecolor(palette[sp]); patch.set_alpha(0.55)
    for i, sp in enumerate(disp_order, start=1):
        y = pg[name][sp_of == sp].values
        x = np.random.normal(i, 0.06, size=len(y))
        axA.scatter(x, y, s=5, color=palette[sp], alpha=0.55, edgecolors="none", zorder=3)
    axA.set_xticks(range(1, len(disp_order) + 1))
    axA.set_xticklabels(disp_order, rotation=40, ha="right", style="italic", fontsize=8)
    axA.set_ylim(0, ymax * 1.05)
    axA.grid(axis="y", ls=":", alpha=0.4)
    axA.set_title(f"{name}\n{subtitle}", fontweight="bold", fontsize=11.5, pad=14)
    axA.text(-0.02, 1.10, f"({letters[0][col]})", transform=axA.transAxes,
             fontweight="bold", fontsize=12, ha="right", va="bottom")
    if col == 0:
        axA.set_ylabel("Defence-system genes\nper genome")

    # ---- row 2: heatmap of top consensus subtypes (0 -> white) ----
    axB = axes[1, col]
    disp = ma.masked_where(hm[name] == 0, np.sqrt(hm[name]))  # mask true absences
    im = axB.imshow(disp, aspect="auto", cmap=CMAP, vmin=0, vmax=sqrt_vmax)
    axB.set_xticks(range(len(disp_order)))
    axB.set_xticklabels(disp_order, rotation=40, ha="right", style="italic", fontsize=8)
    axB.set_yticks(range(len(top_subtypes)))
    axB.set_yticklabels(top_subtypes if col == 0 else [], fontsize=8)
    # light grid so white (absent) cells remain visible as part of the matrix
    axB.set_xticks(np.arange(-.5, len(disp_order), 1), minor=True)
    axB.set_yticks(np.arange(-.5, len(top_subtypes), 1), minor=True)
    axB.grid(which="minor", color="#dddddd", lw=0.5)
    axB.tick_params(which="minor", length=0)
    axB.text(-0.02, 1.015, f"({letters[1][col]})", transform=axB.transAxes,
             fontweight="bold", fontsize=12, ha="right", va="bottom")
    thr = sqrt_vmax * 0.55
    for i in range(len(top_subtypes)):
        for j in range(len(disp_order)):
            v = hm[name][i, j]
            if v >= 1.0:  # only label the prominent cells, keeps panels legible
                axB.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=5.5,
                         color="black" if np.sqrt(v) > thr else "white")

cb = fig2.colorbar(im, ax=axes[1, :].tolist(), fraction=0.03, pad=0.015)
cb.set_label(r"$\sqrt{\mathrm{mean\ copies\ per\ genome}}$")
cb.ax.text(0.5, -0.04, "white =\nnot detected\n(No_hit)", transform=cb.ax.transAxes,
           ha="center", va="top", fontsize=7.5, style="italic")
fig2.suptitle("Defence-system landscape: full Ptolemaea consensus vs. BLAST-only, "
              "PADLOC and DefenseFinder alone\n"
              f"({n_genomes} genomes; per-genome distributions (A–D) share a y-axis; "
              f"heatmaps (E–H) of the top {TOP_N} consensus subtypes share a colour scale)",
              fontweight="bold", fontsize=12.5)
for ext in ("pdf", "png"):
    fig2.savefig(os.path.join(ROOT, "figures", f"fig2_landscape.{ext}"),
                 dpi=200, bbox_inches="tight")

# ---------------------------------------------------------------------------
# FIGURE 3 — concordance
#   (A) 3-set Venn of source membership (PADLOC / DefenseFinder / BLAST)
#   (B) status composition by species (stacked proportions)
# ---------------------------------------------------------------------------
fig3, (axV, axS) = plt.subplots(
    1, 2, figsize=(11, 5.2), gridspec_kw={"width_ratios": [1.0, 1.1]})

# Panel A — non-proportional 3-circle Venn
axV.set_aspect("equal"); axV.axis("off")
axV.set_xlim(-1.7, 1.7); axV.set_ylim(-1.7, 1.9)
r = 1.0
centres = {"P": (-0.5, 0.35), "D": (0.5, 0.35), "B": (0.0, -0.5)}
colours = {"P": "#4C72B0", "D": "#DD8452", "B": "#55A868"}
for k, (cx, cy) in centres.items():
    axV.add_patch(Circle((cx, cy), r, alpha=0.35, color=colours[k], ec="none"))
# region label positions
pos = {
    "P_only": (-1.0, 0.75), "D_only": (1.0, 0.75), "B_only": (0.0, -1.15),
    "PD": (0.0, 0.95), "PB": (-0.62, -0.35), "DB": (0.62, -0.35), "PDB": (0.0, 0.05),
}
for k, (x, y) in pos.items():
    axV.text(x, y, f"{regions[k]:,}", ha="center", va="center", fontsize=9,
             fontweight="bold")
axV.text(-0.95, 1.55, "PADLOC", color=colours["P"], fontweight="bold", ha="center")
axV.text(0.95, 1.55, "DefenseFinder", color=colours["D"], fontweight="bold", ha="center")
axV.text(0.0, -1.62, "BLAST (bidirectional)", color=colours["B"], fontweight="bold",
         ha="center")
axV.set_title("A   Source membership of consensus annotations",
              loc="left", fontweight="bold")

# Panel B — status composition by species (proportions, stacked)
status_order = ["AGREE", "RESOLVED", "SINGLE", "BLAST", "MAPPING", "CONFLICT"]
status_cols = {
    "AGREE": "#2c7bb6", "RESOLVED": "#abd9e9", "SINGLE": "#ffffbf",
    "BLAST": "#fdae61", "MAPPING": "#d7191c", "CONFLICT": "#7b3294",
}
props = []
for sp in disp_order:
    sub = ALL[ALL["species"] == sp]["status"].value_counts()
    tot = sub.sum()
    props.append([sub.get(st, 0) / tot * 100 for st in status_order])
props = np.array(props)
ypos = np.arange(len(disp_order))
left = np.zeros(len(disp_order))
for k, st in enumerate(status_order):
    axS.barh(ypos, props[:, k], left=left, color=status_cols[st], label=st,
             edgecolor="white", height=0.7)
    left += props[:, k]
axS.set_yticks(ypos); axS.set_yticklabels(disp_order, style="italic")
axS.invert_yaxis()
axS.set_xlabel("Percentage of consensus annotations")
axS.set_xlim(0, 100)
axS.set_title("B   How each consensus call was reached",
              loc="left", fontweight="bold")
axS.legend(ncol=3, fontsize=7.5, loc="upper center",
           bbox_to_anchor=(0.5, -0.12), frameon=False)
fig3.suptitle("Concordance of defence-system calls between sources",
              fontweight="bold", y=1.01)
fig3.tight_layout()
for ext in ("pdf", "png"):
    fig3.savefig(os.path.join(ROOT, "figures", f"fig3_concordance.{ext}"),
                 dpi=200, bbox_inches="tight")

print("\nWrote: results_stats.txt, figures/fig2_landscape.{pdf,png}, "
      "figures/fig3_concordance.{pdf,png}")
