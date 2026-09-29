#!/usr/bin/env python3
# iTOL connection files for the spacer to prophage links (Fig 5.11-5.13).
# Usage: python3 scripts/make_itol_crispr_linkages.py 2 3 5

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

DATA = Path("data")
HITS = DATA / "Protospacer_hits.csv"
COMM = DATA / "genome_community_mapping_49_BAC.csv"
OUTDIR = Path("figures/itol_linkages")

GREY, RED, BLUE = "#cccccc", "#e31a1c", "#1f78b4"

HEADER = """DATASET_CONNECTION
SEPARATOR COMMA

DATASET_LABEL,CRISPR spacer -> protospacer (community {comm})
COLOR,#ff0000

DRAW_ARROWS,1
ARROW_SIZE,20
MAXIMUM_LINE_WIDTH,10
CURVE_ANGLE,0
CENTER_CURVES,1
ALIGN_TO_LABELS,0

LEGEND_TITLE,Spacer targeting
LEGEND_SHAPES,1,1,1
LEGEND_COLORS,{red},{blue},{grey}
LEGEND_LABELS,Inside community targets outside,Outside targets inside community,Within community

DATA
"""


def load_pairs(hits_path):
    """Count spacer->prophage hits per (source genome, target genome) pair."""
    pairs = Counter()
    with open(hits_path, newline="") as fh:
        rd = csv.reader(fh)
        next(rd, None)
        for row in rd:
            if len(row) < 2:
                continue
            src = row[0].split("|")[0].strip()
            tgt = row[1].split("_")[0].strip()
            if src and tgt:
                pairs[(src, tgt)] += 1
    return pairs


def load_communities(comm_path):
    """Genome ID -> bacterial community. Strips the .fna suffix."""
    g2c = {}
    with open(comm_path, newline="") as fh:
        for row in csv.DictReader(fh):
            gid = row["Genome_ID"].strip()
            if gid.endswith(".fna"):
                gid = gid[:-4]
            g2c[gid] = row["Community"].strip()
    return g2c


def write_linkage(comm, pairs, g2c, outdir):
    members = {g for g, c in g2c.items() if c == str(comm)}
    if not members:
        return None, 0

    rows = []
    for (src, tgt), n in pairs.items():
        s_in, t_in = src in members, tgt in members
        if not (s_in or t_in):
            continue
        colour = GREY if (s_in and t_in) else (RED if s_in else BLUE)
        rows.append(f"{src},{tgt},{1.0 + 0.45 * n:.2f},{colour},normal")

    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"Linkage_community_{comm}.txt"
    with open(out, "w") as fh:
        fh.write(HEADER.format(comm=comm, red=RED, blue=BLUE, grey=GREY))
        fh.write("\n".join(rows) + "\n")
    return out, len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("communities", nargs="*", help="community IDs, e.g. 3 5 2")
    ap.add_argument("--all", action="store_true", help="every community with links")
    args = ap.parse_args()

    pairs = load_pairs(HITS)
    g2c = load_communities(COMM)
    print(f"{sum(pairs.values()):,} hits across {len(pairs):,} genome pairs")
    print(f"{len(g2c):,} genomes in {len(set(g2c.values())):,} communities")

    targets = sorted(set(g2c.values()), key=lambda x: int(x)) if args.all \
        else (args.communities or ["3", "5", "2"])

    for comm in targets:
        out, n = write_linkage(comm, pairs, g2c, OUTDIR)
        if out and n:
            print(f"  community {comm:>4}: {n:>6,} connections -> {out}")


if __name__ == "__main__":
    main()
