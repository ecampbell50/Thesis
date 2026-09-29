#!/usr/bin/env python3
# CRISPR targeting from the spacer vs provirus BLAST hits (hpc/07).
# Writes tables/crispr_*.csv.

import os, sys
import pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
RES  = "data"
COE  = "data"
HITS = os.path.join("data", "spacer_vs_provirus.tsv")

if not os.path.exists(HITS):
    sys.exit(f"Not found: {HITS}\nRun crispr_spacer_pipeline.sh on the HPC and scp the tsv here first.")

cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
        "qlen","qstart","qend","sstart","send","evalue","bitscore"]
hits = pd.read_csv(HITS, sep="\t", names=cols, dtype=str)
for c in ["pident","length","mismatch","qlen"]:
    hits[c] = pd.to_numeric(hits[c])

# spacer host genome = text before the first '|' in the spacer id
hits["spacer_host"] = hits["qseqid"].str.split("|").str[0]
hits = hits.rename(columns={"sseqid": "prophage_id"})

# provirus -> (pro_comm, origin host, origin bac_comm)
long = pd.read_csv(os.path.join("data", "prophage_host_community_long.csv"),
                   dtype=str).rename(columns={"host": "origin_host", "bac_comm": "origin_baccomm"})
hits = hits.merge(long[["prophage_id","pro_comm","origin_host","origin_baccomm"]],
                  on="prophage_id", how="left")

# spacer host -> its BAC community
bac = pd.read_csv(os.path.join("data", "genome_community_mapping_49_BAC.csv"), dtype=str)
bac["host"] = bac["Genome_ID"].str.replace(".fna", "", regex=False)
bac["spacer_host_baccomm"] = "BAC_com_" + bac["Community"].astype(str)
hits = hits.merge(bac[["host","spacer_host_baccomm"]],
                  left_on="spacer_host", right_on="host", how="left").drop(columns=["host"])

# classify each targeting event
def kind(r):
    if r["spacer_host"] == r["origin_host"]:           return "self"
    if r["spacer_host_baccomm"] == r["origin_baccomm"]: return "within_lineage"
    return "cross_lineage"
hits["targeting_type"] = hits.apply(kind, axis=1)
hits.to_csv(os.path.join("tables", "crispr_targeting_events.csv"), index=False)

# per prophage-community immunity pressure
pressure = (hits.groupby("pro_comm")
            .agg(n_spacer_hits=("qseqid","size"),
                 n_targeting_hosts=("spacer_host","nunique"),
                 n_targeting_lineages=("spacer_host_baccomm","nunique"))
            .sort_values("n_targeting_lineages", ascending=False))
pressure.to_csv(os.path.join("tables", "crispr_procomm_pressure.csv"))

# immunity vs carriage, done PROPERLY: enumerate the full universe of
# (spacer-bearing lineage  x  real prophage community) pairs so the true-negative
# cell exists, then Fisher. Without the negatives the 2x2 is a base-rate mirage:
# carriage is rare overall, so "targeted-but-not-carried" is huge regardless.
from scipy.stats import fisher_exact
real = list(pd.read_csv(os.path.join("data", "host_by_procomm.csv"),
                        index_col=0, nrows=0).columns)            # the 85 real communities
evr  = hits[hits.pro_comm.isin(real)]
lor  = long[long.pro_comm.isin(real)]
L = sorted(evr.spacer_host_baccomm.dropna().unique())            # lineages where immunity is observable
T = set(map(tuple, evr[["spacer_host_baccomm","pro_comm"]].dropna().drop_duplicates().values))
C = set(map(tuple, lor[["origin_baccomm","pro_comm"]].dropna().drop_duplicates().values))
a = b = c = d = 0
for lin in L:
    for pc in real:
        t = (lin, pc) in T; k = (lin, pc) in C
        if   t and k:        a += 1
        elif t and not k:    b += 1
        elif (not t) and k:  c += 1
        else:                d += 1
OR, p = fisher_exact([[a, b], [c, d]])
tab = pd.DataFrame([[a, b], [c, d]],
                   index=["targets=1", "targets=0"], columns=["carries=1", "carries=0"])
tab.to_csv(os.path.join("tables", "crispr_immunity_vs_carriage.csv"))
with open(os.path.join(COE, "crispr_immunity_vs_carriage_stats.txt"), "w") as fh:
    fh.write(f"universe: {len(L)} spacer-bearing lineages x {len(real)} real prophage communities\n")
    fh.write(f"P(carry|target)={a/(a+b):.3f}  P(carry|no target)={c/(c+d):.3f}\n")
    fh.write(f"Fisher OR={OR:.3f}  p={p:.3e}  (OR>1 => spacer presence POSITIVELY associated with carriage)\n")

print(f"events: {len(hits)}  | self={ (hits.targeting_type=='self').sum() }"
      f"  within={ (hits.targeting_type=='within_lineage').sum() }"
      f"  cross={ (hits.targeting_type=='cross_lineage').sum() }")
print("\nTop targeted prophage communities:")
print(pressure.head(10))
print(f"\nimmunity vs carriage (proper 2x2, {len(L)} lineages x {len(real)} communities):")
print(tab.to_string())
print(f"  P(carry|target)={a/(a+b):.1%}  P(carry|no target)={c/(c+d):.1%}  Fisher OR={OR:.2f}  p={p:.2e}")
