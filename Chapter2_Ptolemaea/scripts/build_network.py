#!/usr/bin/env python3
# ======================================================================
# Thesis Chapter 2, Figure 2.1
# Defence-system co-annotation network panel (run with --figure to draw it).
# Original location: github.com/ecampbell50/Chp4_Ptolemaea scripts/build_network.py
# ======================================================================
"""
build_network.py — cross-tool co-annotation / conflict network for Ptolemaea.

For every predicted protein, PADLOC, DefenseFinder and the forward/reverse BLAST
searches may each emit a (raw) system name. When two or more *different* names are
applied to the same protein, that is an inter-tool disagreement. Aggregated across
all 700 genomes, these co-annotations form a network:
    * node  = a raw tool-level system label
    * edge  = two labels applied to the same protein by different sources
    * weight= number of proteins showing that co-annotation

Usage:
    python3 scripts/build_network.py            # exploration printout
    python3 scripts/build_network.py --figure   # also render figures/fig4_network
"""
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPECIES = ["Ecoli", "Efaecium", "Saureus", "Kpneumoniae",
           "Abaumanii", "Paeruginosa", "Enterobacter"]
NULL = "No_hit"
SRC_COLS = ["padloc_original", "deffind_original", "fwd_blast", "rev_blast"]


def clean(s):
    return s.replace("CRISP-Cas", "CRISPR-Cas") if isinstance(s, str) else s


def load():
    frames = []
    for tag in SPECIES:
        a = pd.read_csv(os.path.join(ROOT, "data", f"{tag}_annotations.csv"), dtype=str).fillna(NULL)
        for c in SRC_COLS:
            a[c] = a[c].map(clean)
        frames.append(a)
    return pd.concat(frames, ignore_index=True)


CONTESTED = {"RESOLVED", "CONFLICT", "MAPPING"}  # statuses that mark genuine disagreement


def build(df, exclude_agree=False):
    """Return (node_freq, edge_weight, label_source, edge_contested)."""
    node_freq = Counter()
    edge_w = Counter()
    edge_contested = Counter()              # edge -> # proteins with a contested status
    src_col_counts = defaultdict(Counter)   # label -> Counter({tool: n})
    colmap = {"padloc_original": "PADLOC", "deffind_original": "DefenseFinder",
              "fwd_blast": "BLAST", "rev_blast": "BLAST"}

    sub = df[df["status"] != "AGREE"] if exclude_agree else df
    for _, row in sub[["status"] + SRC_COLS].iterrows():
        labels = {}
        for c in SRC_COLS:
            v = row[c]
            if v != NULL:
                src_col_counts[v][colmap[c]] += 1
                labels.setdefault(v, set()).add(colmap[c])
        uniq = list(labels.keys())
        for lab in uniq:
            node_freq[lab] += 1
        if len(uniq) >= 2:
            contested = row["status"] in CONTESTED
            for a, b in combinations(sorted(uniq), 2):
                edge_w[(a, b)] += 1
                if contested:
                    edge_contested[(a, b)] += 1

    label_source = {}
    for lab, cnt in src_col_counts.items():
        tools = set(cnt.keys())
        label_source[lab] = next(iter(tools)) if len(tools) == 1 else "shared"
    return node_freq, edge_w, label_source, edge_contested


TOOL_COL = {"PADLOC": "#4C72B0", "DefenseFinder": "#DD8452",
            "BLAST": "#55A868", "shared": "#9b9b9b"}


def packed_layout(G, nx, np, pad=1.0):
    """Lay out each connected component locally, then shelf-pack the tiles so
    nothing overlaps and there is no wasted whitespace."""
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    tiles = []  # (nodes, local_pos, side)
    for comp in comps:
        sub = G.subgraph(comp)
        n = len(comp)
        if n == 2:
            a, b = list(comp)
            lp = {a: np.array([0.0, 0.0]), b: np.array([1.0, 0.0])}
        else:
            lp = nx.spring_layout(sub, weight="weight", k=1.3,
                                  iterations=500, seed=42)
        xs = [p[0] for p in lp.values()]; ys = [p[1] for p in lp.values()]
        minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)
        sx, sy = (maxx - minx) or 1.0, (maxy - miny) or 1.0
        side = 1.7 * np.sqrt(n)
        norm = {nd: np.array([(p[0] - minx) / sx, (p[1] - miny) / sy]) * side
                for nd, p in lp.items()}
        tiles.append((norm, side))

    # shelf packing: order by side desc, wrap rows at a target width
    target_w = max(t[1] for t in tiles) * 0  # placeholder
    total = sum(t[1] + pad for t in tiles)
    target_w = max(total / 4.5, max(t[1] for t in tiles))
    pos = {}
    cx, cy, row_h = 0.0, 0.0, 0.0
    for norm, side in tiles:
        if cx > 0 and cx + side > target_w:
            cx = 0.0; cy -= (row_h + pad); row_h = 0.0
        for nd, p in norm.items():
            pos[nd] = (cx + p[0], cy + p[1])
        cx += side + pad
        row_h = max(row_h, side)
    return pos


def render(node_freq, edge_w, label_source, edge_contested, thresh, out):
    import numpy as np
    import networkx as nx
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    edges = [(a, b, w) for (a, b), w in edge_w.items() if w >= thresh]
    G = nx.Graph()
    for a, b, w in edges:
        contested = edge_contested[(a, b)] / w >= 0.5
        G.add_edge(a, b, weight=w, contested=contested)

    # focus the figure: keep clusters of >=3 tool-names, plus any pair whose
    # co-annotation was contested (a genuine conflict); drop pure synonym pairs.
    keep = set()
    for comp in nx.connected_components(G):
        sub = G.subgraph(comp)
        if len(comp) >= 3 or any(d["contested"] for *_, d in sub.edges(data=True)):
            keep |= set(comp)
    G = G.subgraph(keep).copy()
    print(f"Rendering network: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges (weight >= {thresh}; "
          "pure synonym pairs omitted)")

    pos = packed_layout(G, nx, np)

    fig, ax = plt.subplots(figsize=(15, 12))
    ax.axis("off")

    grey = [(u, v) for u, v, d in G.edges(data=True) if not d["contested"]]
    red = [(u, v) for u, v, d in G.edges(data=True) if d["contested"]]
    wmax = max(d["weight"] for *_, d in G.edges(data=True))

    def widths(elist):
        return [0.5 + 5.0 * (G[u][v]["weight"] / wmax) for u, v in elist]

    nx.draw_networkx_edges(G, pos, edgelist=grey, width=widths(grey),
                           edge_color="#c9c9c9", alpha=0.8, ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=red, width=widths(red),
                           edge_color="#c0392b", alpha=0.9, ax=ax)

    nodes = list(G.nodes())
    sizes = [55 + 16 * np.sqrt(node_freq[n]) for n in nodes]
    colors = [TOOL_COL.get(label_source.get(n, "shared"), "#9b9b9b") for n in nodes]
    nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_size=sizes,
                           node_color=colors, edgecolors="white",
                           linewidths=0.8, ax=ax)

    # label every node, nudged just above its marker to reduce overlap
    label_pos = {n: (x, y + 0.16) for n, (x, y) in pos.items()}
    nx.draw_networkx_labels(G, label_pos, labels={n: n for n in nodes},
                            font_size=6.3, ax=ax,
                            bbox=dict(boxstyle="round,pad=0.1", fc="white",
                                      ec="none", alpha=0.7))

    legend = [
        Line2D([0], [0], marker="o", color="w", label="PADLOC",
               markerfacecolor=TOOL_COL["PADLOC"], markersize=11),
        Line2D([0], [0], marker="o", color="w", label="DefenseFinder",
               markerfacecolor=TOOL_COL["DefenseFinder"], markersize=11),
        Line2D([0], [0], marker="o", color="w", label="BLAST",
               markerfacecolor=TOOL_COL["BLAST"], markersize=11),
        Line2D([0], [0], marker="o", color="w", label="shared vocabulary",
               markerfacecolor=TOOL_COL["shared"], markersize=11),
        Line2D([0], [0], color="#c0392b", lw=3,
               label="contested (RESOLVED/CONFLICT/MAPPING)"),
        Line2D([0], [0], color="#cfcfcf", lw=3, label="concordant (naming variant)"),
    ]
    ax.legend(handles=legend, loc="lower left", fontsize=9, frameon=True,
              framealpha=0.9, title="node = tool of origin   ·   edge = co-annotation")
    ax.set_title("Cross-tool co-annotation network of defence-system calls\n"
                 f"(700 genomes; edges shown for co-annotations on ≥ {thresh} proteins; "
                 "node size ∝ frequency, edge width ∝ co-occurrence)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=200, bbox_inches="tight")

    # contested edges that survive the filter -> for the manuscript text
    print("\nContested (red) edges in the figure, by weight:")
    for (a, b), w in sorted(edge_w.items(), key=lambda kv: -kv[1]):
        if w >= thresh and edge_contested[(a, b)] / w >= 0.5:
            print(f"   {w:>4}  {a} [{label_source.get(a)}]  <->  "
                  f"{b} [{label_source.get(b)}]   "
                  f"(contested {edge_contested[(a,b)]}/{w})")


def report(node_freq, edge_w, label_source, title):
    print("=" * 70); print(title); print("=" * 70)
    print(f"nodes: {len(node_freq):,}   edges: {len(edge_w):,}")
    print("\nTop 30 co-annotation edges (weight = # proteins):")
    for (a, b), w in edge_w.most_common(30):
        print(f"   {w:>5}  {a} [{label_source.get(a,'?')}]  <->  {b} [{label_source.get(b,'?')}]")
    print()


if __name__ == "__main__":
    df = load()
    print(f"Loaded {len(df):,} protein annotations across {df['genome_id'].nunique()} genomes\n")
    nf, ew, ls, ec = build(df, exclude_agree=False)
    if "--figure" in sys.argv:
        thresh = 40
        for a in sys.argv:
            if a.startswith("--thresh="):
                thresh = int(a.split("=")[1])
        render(nf, ew, ls, ec, thresh, os.path.join(ROOT, "figures", "fig4_network"))
    else:
        report(nf, ew, ls, "ALL proteins (incl. AGREE)")
