#!/usr/bin/env python3
# Per-community summary tables for BipartiteNetwork_Analysis.ipynb.

import os, re
import numpy as np, pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

RES = "Results"
OUT = os.path.join(RES, "community_questions")
os.makedirs(OUT, exist_ok=True)

# ---- community maps (string-safe; keep .fna on IDs for edge matching) -------
bac = pd.read_csv(os.path.join(RES, "data/genome_community_mapping_49_BAC.csv"), dtype=str)
pro = pd.read_csv(os.path.join(RES, "data/genome_community_mapping_44_PRO.csv"), dtype=str)
bac.columns = ["g", "c"]; pro.columns = ["g", "c"]
bac["bac_comm"] = "BAC_com_" + bac["c"]
pro["pro_comm"] = "PRO_com_" + pro["c"]
bac_set, pro_set = set(bac["g"]), set(pro["g"])
g2bac = dict(zip(bac["g"], bac["bac_comm"]))
g2pro = dict(zip(pro["g"], pro["pro_comm"]))

# community sizes (genomes per community)
bac_size = bac.groupby("bac_comm").size().rename("n_genomes")
pro_size = pro.groupby("pro_comm").size().rename("n_prophages")

# ---- host edges (BAC<->PRO) from the bipartite edgetable --------------------
edf = pd.read_csv(os.path.join(RES, "data/Bipartite_BAC_PRO_genome_edgetable.csv"), dtype=str)
edf.columns = ["Source", "Target", "Value"]
host = edf[(edf.Source.isin(bac_set) & edf.Target.isin(pro_set)) |
           (edf.Source.isin(pro_set) & edf.Target.isin(bac_set))].copy()
# normalise so col 'bg' = bacterial genome, 'pg' = prophage genome
host["bg"] = np.where(host.Source.isin(bac_set), host.Source, host.Target)
host["pg"] = np.where(host.Source.isin(pro_set), host.Source, host.Target)
host["bac_comm"] = host["bg"].map(g2bac)
host["pro_comm"] = host["pg"].map(g2pro)
pairs = host[["bac_comm", "pro_comm"]].drop_duplicates()   # community-community links

# ============================================================================
# Q1  prophage community size vs # bacterial-community connections
# ============================================================================
q1 = (pairs.groupby("pro_comm").size().rename("n_baccomm_connections")
      .to_frame().join(pro_size, how="right").fillna({"n_baccomm_connections": 0}))
q1["n_baccomm_connections"] = q1["n_baccomm_connections"].astype(int)
q1 = q1.reset_index()[["pro_comm", "n_prophages", "n_baccomm_connections"]]
q1.to_csv(os.path.join(OUT, "data/q1_procomm_promiscuity.csv"), index=False)

# ============================================================================
# Q2  prophage community: mean length vs size, coloured by SD of length
#     (length parsed from the provirus coords in each prophage ID)
# ============================================================================
def plen(pid):
    m = re.search(r"provirus-(\d+)-(\d+)", pid)
    return abs(int(m.group(2)) - int(m.group(1))) if m else np.nan
pro["length"] = pro["g"].map(plen)
q2 = (pro.groupby("pro_comm")["length"]
        .agg(n_prophages="size", mean_length="mean", sd_length="std", median_length="median")
        .reset_index())
q2.to_csv(os.path.join(OUT, "data/q2_procomm_lengths.csv"), index=False)

# ============================================================================
# Q3  bacterial community size vs # prophage communities (and per-member rate)
# ============================================================================
q3 = (pairs.groupby("bac_comm").size().rename("n_procomm_connections")
      .to_frame().join(bac_size, how="right").fillna({"n_procomm_connections": 0}))
q3["n_procomm_connections"] = q3["n_procomm_connections"].astype(int)
q3["procomm_per_member"] = q3["n_procomm_connections"] / q3["n_genomes"]
q3 = q3.reset_index()[["bac_comm", "n_genomes", "n_procomm_connections", "procomm_per_member"]]
q3.to_csv(os.path.join(OUT, "data/q3_baccomm_procoms.csv"), index=False)

# ============================================================================
# Q4  bacterial-community defence profiles (binomial exact, BAC-comm x subtype)
#     mirrors the Chapter-1 enrichment table so it can feed the same heatmap.
# ============================================================================
d = pd.read_csv(os.path.join(RES, "data/Ssuis_DefenceMatrix_CLEAN_FINAL.csv"), dtype=str)
d = d.rename(columns={d.columns[0]: "genome_id"}).set_index("genome_id")
d = (d.apply(pd.to_numeric, errors="coerce").fillna(0) > 0).astype(int)
# collapse columns to SUBTYPE level
def subt(c):
    m = re.search(r"SUBTYPE#([^#]+)#", c); return m.group(1) if m else c
d = d.groupby(subt, axis=1).max()
# join genome -> bac community (strip .fna from mapping IDs to match defence IDs)
gid2bac = {g.replace(".fna", ""): c for g, c in g2bac.items()}
d = d.loc[d.index.intersection(gid2bac.keys())]
comm_of = pd.Series({g: gid2bac[g] for g in d.index})
pop_freq = d.mean(axis=0)                       # overall carriage per subtype
rows = []
for comm, idx in comm_of.groupby(comm_of).groups.items():
    sub = d.loc[idx]; n = len(idx)
    for st in d.columns:
        k = int(sub[st].sum()); p = float(pop_freq[st])
        if p <= 0 or p >= 1:                     # skip universal/absent subtypes
            continue
        exp = n * p
        bp = binomtest(k, n, p, alternative="two-sided").pvalue
        rows.append(dict(bac_comm=comm, defence_subtype=st, n_carrying=k, n_in_comm=n,
                         pop_freq=p, expected=exp,
                         log2OE=np.log2((k + 0.5) / (exp + 0.5)), binom_p=bp))
q4 = pd.DataFrame(rows)
q4["q_BH"] = multipletests(q4["binom_p"], method="fdr_bh")[1]
q4["significant"] = q4["q_BH"] < 0.05
q4.to_csv(os.path.join(OUT, "data/q4_baccomm_defence_enrichment.csv"), index=False)

print("wrote to", OUT)
for f, df in [("q1_procomm_promiscuity", q1), ("q2_procomm_lengths", q2),
              ("q3_baccomm_procoms", q3), ("q4_baccomm_defence_enrichment", q4)]:
    print(f"  {f}.csv  ({len(df)} rows)")
print(f"\nQ4: {int(q4['significant'].sum())} significant bac-comm x subtype enrichments (q<0.05)")
