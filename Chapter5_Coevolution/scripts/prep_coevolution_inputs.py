#!/usr/bin/env python3
"""
prep_coevolution_inputs.py
==========================
Builds the clean input matrices for the two local coevolution analyses:

  Move 1  Network structure (nestedness / modularity)   -> analyse_network_structure.R
  Move 2  Pagel correlated-evolution test               -> analyse_pagel.R

All outputs land in Results/coevolution/.  Wrangling lives here (pandas);
the statistics live in the R scripts.  genome_id is ALWAYS read as text
(the 1307.1890 -> 1307.189 float trap).

Inputs
------
  Results/Ssuis_DefenceMatrix_CLEAN_FINAL.csv             genome_id x 149 system COUNT columns
  Results/prophage_defence/prophage_host_community_long.csv   prophage_id, pro_comm, host, bac_comm
  Results/prophage_defence/community_defence_carriage.csv     pro_comm, system, frac_hosts_carrying, flag, ...
  Results/prophage_defence/community_defence_enrichment.csv   pro_comm, system, direction, q_value, significant
  iToL/SsuisPhylo_FastTree_11Sep24.treefile              (tip labels only, read here to report overlap)

Outputs (Results/coevolution/)
------------------------------
  host_by_procomm.csv          hosts x prophage-communities, binary  (infection network)
  procomm_by_defence_frac.csv  prophage-communities x defence systems, weighted by frac carriage
  host_by_defence_binary.csv   hosts x defence systems, presence/absence
  pagel_pairs.csv              shortlist of (pro_comm, defence_system) pairs to test phylogenetically
  prep_summary.txt             diagnostics (sizes, overlaps, thresholds)
"""

import os
import re
import sys
import pandas as pd
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # Chapter5_Coevolution/
RES  = ROOT
PD_DIR = ROOT
OUT  = ROOT
os.makedirs(OUT, exist_ok=True)

# ---- tunable thresholds -----------------------------------------------------
MIN_COMM_HOSTS = 5     # only treat prophage communities with >= this many distinct hosts as "real"
                       # (matches min_summary_sz in plot_bipartite_network.R)
N_PAGEL_PAIRS  = 20    # how many top enriched + top depleted pairs to queue for Pagel
# -----------------------------------------------------------------------------

log_lines = []
def log(msg=""):
    print(msg)
    log_lines.append(str(msg))

# ---------------------------------------------------------------------------
# 1. Defence matrix -> system-level presence/absence per host genome
# ---------------------------------------------------------------------------
log("== Defence matrix ==")
dm = pd.read_csv(os.path.join(RES, "data/Ssuis_DefenceMatrix_CLEAN_FINAL.csv"),
                 dtype={"genome_id": str})
dm = dm.set_index("genome_id")

# Column names look like  TYPE#AVAST#SUBTYPE#AVAST_II#OUTCOME#Abi
# Collapse to system label  "TYPE:SUBTYPE"  to match the carriage/enrichment 'system' field.
def col_to_system(col):
    m = re.match(r"TYPE#(.*?)#SUBTYPE#(.*?)#OUTCOME#(.*)", col)
    if not m:
        return col  # leave anything unexpected untouched
    typ, sub = m.group(1), m.group(2)
    return f"{typ}:{sub}"

sys_of_col = {c: col_to_system(c) for c in dm.columns}
# presence = count > 0, then collapse columns that map to the same system (any present)
pres = (dm > 0).astype(int)
pres.columns = [sys_of_col[c] for c in pres.columns]
host_by_defence = pres.groupby(level=0, axis=1).max()   # system-level presence/absence
log(f"  genomes: {host_by_defence.shape[0]}   raw count cols: {dm.shape[1]}   "
    f"-> systems: {host_by_defence.shape[1]}")

# ---------------------------------------------------------------------------
# 2. prophage_host_community_long -> host x prophage-community incidence
# ---------------------------------------------------------------------------
log("\n== Prophage host-community links ==")
long = pd.read_csv(os.path.join(PD_DIR, "data/prophage_host_community_long.csv"),
                   dtype={"host": str, "pro_comm": str, "bac_comm": str})
log(f"  links: {len(long)}   distinct hosts w/ prophage: {long['host'].nunique()}   "
    f"distinct pro_comms: {long['pro_comm'].nunique()}")

# size of each prophage community = number of DISTINCT host genomes it infects
comm_hostcount = long.groupby("pro_comm")["host"].nunique().sort_values(ascending=False)
real_comms = comm_hostcount[comm_hostcount >= MIN_COMM_HOSTS].index.tolist()
log(f"  prophage communities with >= {MIN_COMM_HOSTS} distinct hosts: {len(real_comms)} "
    f"(of {len(comm_hostcount)} total)")

# host x pro_comm incidence (restricted to real communities)
sub = long[long["pro_comm"].isin(real_comms)]
host_by_procomm = (pd.crosstab(sub["host"], sub["pro_comm"]) > 0).astype(int)
log(f"  host_by_procomm incidence: {host_by_procomm.shape[0]} hosts x "
    f"{host_by_procomm.shape[1]} prophage communities")

# lineage-level version: collapse hosts within a BAC community (pseudoreplication control).
# A BAC community 'carries' a prophage community if ANY member host does.
baccomm_by_procomm = (pd.crosstab(sub["bac_comm"], sub["pro_comm"]) > 0).astype(int)
log(f"  baccomm_by_procomm incidence: {baccomm_by_procomm.shape[0]} BAC communities x "
    f"{baccomm_by_procomm.shape[1]} prophage communities")

# ---------------------------------------------------------------------------
# 3. community_defence_carriage -> pro_comm x defence (weighted by frac carriage)
# ---------------------------------------------------------------------------
log("\n== Prophage-community defence carriage ==")
carr = pd.read_csv(os.path.join(PD_DIR, "data/community_defence_carriage.csv"))
carr = carr[carr["pro_comm"].isin(real_comms)]
procomm_by_defence = carr.pivot_table(index="pro_comm", columns="system",
                                      values="frac_hosts_carrying", fill_value=0.0)
# drop systems that never appear in any real community and that are universal (no signal)
nz = (procomm_by_defence > 0).sum(axis=0)
keep = nz[(nz > 0) & (nz < procomm_by_defence.shape[0])].index
procomm_by_defence = procomm_by_defence[keep]
log(f"  procomm_by_defence (frac-weighted): {procomm_by_defence.shape[0]} communities x "
    f"{procomm_by_defence.shape[1]} systems (after dropping empty/universal)")

# ---------------------------------------------------------------------------
# 4. Pagel pair shortlist: strongest enriched & depleted associations
# ---------------------------------------------------------------------------
log("\n== Pagel pair shortlist ==")
enr = pd.read_csv(os.path.join(PD_DIR, "data/community_defence_enrichment.csv"))
enr = enr[(enr["significant"] == True) & (enr["pro_comm"].isin(real_comms))].copy()
enr["q_value"] = pd.to_numeric(enr["q_value"], errors="coerce")
pairs = []
for direction in ("enriched", "depleted"):
    d = enr[enr["direction"] == direction].sort_values("q_value").head(N_PAGEL_PAIRS)
    pairs.append(d[["pro_comm", "system", "direction", "q_value",
                    "frac_hosts_carrying", "pop_freq"]])
pagel_pairs = pd.concat(pairs, ignore_index=True)
# only keep pairs whose defence system actually exists as a host-level column
pagel_pairs = pagel_pairs[pagel_pairs["system"].isin(host_by_defence.columns)]
log(f"  queued {len(pagel_pairs)} pairs "
    f"({(pagel_pairs.direction=='enriched').sum()} enriched / "
    f"{(pagel_pairs.direction=='depleted').sum()} depleted)")

# ---------------------------------------------------------------------------
# 5. Tree overlap report (don't prune here; R does that against the actual tree)
# ---------------------------------------------------------------------------
log("\n== Tree overlap ==")
tips = set()
try:
    with open(os.path.join(ROOT, "iToL", "data/SsuisPhylo_FastTree_11Sep24.treefile")) as fh:
        nwk = fh.read()
    tips = set(re.findall(r"[\(,]([^(),:]+):", nwk))
    log(f"  tree tips parsed: {len(tips)}")
    log(f"  tips ∩ defence-matrix genomes: {len(tips & set(host_by_defence.index))}")
    log(f"  tips ∩ hosts-with-prophage:    {len(tips & set(host_by_procomm.index))}")
except Exception as e:
    log(f"  (could not parse tree here: {e}; R will handle the join)")

# ---------------------------------------------------------------------------
# 6. Write outputs
# ---------------------------------------------------------------------------
host_by_procomm.to_csv(os.path.join(OUT, "data/host_by_procomm.csv"))
baccomm_by_procomm.to_csv(os.path.join(OUT, "data/baccomm_by_procomm.csv"))
procomm_by_defence.to_csv(os.path.join(OUT, "data/procomm_by_defence_frac.csv"))
host_by_defence.to_csv(os.path.join(OUT, "data/host_by_defence_binary.csv"))
pagel_pairs.to_csv(os.path.join(OUT, "data/pagel_pairs.csv"), index=False)

with open(os.path.join(OUT, "prep_summary.txt"), "w") as fh:
    fh.write("\n".join(log_lines) + "\n")

log(f"\nWrote 4 matrices + prep_summary.txt to {OUT}")
