#!/usr/bin/env python3
# =============================================================================
# Prophage community  x  host defence-system analysis
#
# Hypothesis (chapter aim):
#   Prophages cluster into communities because they are genomically similar
#   (sourmash). A single prophage community can sit inside hosts that fall in
#   several *different* bacterial communities (the prophages are similar, but
#   their hosts are not similar enough to co-cluster). So: when one prophage
#   community spans multiple bacterial communities, is there something common
#   among those hosts -- e.g. the defence systems they carry?
#
# This script answers the defence-system part:
#   For every prophage community, take the set of distinct HOST genomes of its
#   member prophages and ask which defence systems are ALWAYS present, NEVER
#   present, or statistically enriched / depleted relative to the whole
#   S. suis population (Fisher's exact test, BH-FDR corrected).
#
# Join chain (validated 100% coverage):
#   prophage  --(ID prefix before first '_')-->  host genome
#   host genome --> its bacterial community (genome_community_mapping_49_BAC)
#   host genome --> its defence profile      (Ssuis_DefenceMatrix_CLEAN_FINAL)
#
# IMPORTANT: genome_id is read as TEXT everywhere. As a float, '1307.1890'
# collapses onto '1307.189' (trailing zero dropped) and two distinct BV-BRC
# genomes merge. dtype=str on every load + every merge key keeps them apart.
# =============================================================================

import os
import re
import numpy as np
import pandas as pd
from scipy.stats import binom

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
IN_DIR    = "data"
OUT_DIR   = "data"

DEF_FILE  = os.path.join(IN_DIR, "Ssuis_DefenceMatrix_CLEAN_FINAL.csv")
PRO_FILE  = os.path.join(IN_DIR, "genome_community_mapping_44_PRO.csv")
BAC_FILE  = os.path.join(IN_DIR, "genome_community_mapping_49_BAC.csv")

GRAIN     = "subtype"   # "type"  -> one column per defence system TYPE (CBASS, RM, Thoeris, ...)
                     # "subtype" -> finer: one column per TYPE:SUBTYPE
MIN_HOSTS = 3        # only run enrichment tests for communities with >= this many distinct hosts
FDR_ALPHA = 0.05     # BH-FDR significance threshold flagged in the output

os.makedirs(OUT_DIR, exist_ok=True)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def benjamini_hochberg(pvals):
    """Return BH-FDR q-values for a 1D array of p-values (NaNs preserved)."""
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    qv = ranked * n / (np.arange(1, n + 1))
    qv = np.minimum.accumulate(qv[::-1])[::-1]   # enforce monotonicity
    out = np.empty(n)
    out[order] = np.clip(qv, 0, 1)
    q[ok] = out
    return q


def parse_def_columns(cols):
    """TYPE#x#SUBTYPE#y#OUTCOME#z -> dict(type, subtype, outcome)."""
    meta = {}
    for c in cols:
        m = re.match(r"TYPE#(.*)#SUBTYPE#(.*)#OUTCOME#(.*)$", c)
        if m:
            meta[c] = dict(type=m.group(1), subtype=m.group(2), outcome=m.group(3))
        else:
            meta[c] = dict(type=c, subtype=c, outcome="NA")
    return meta


# ----------------------------------------------------------------------------
# 1. Load
# ----------------------------------------------------------------------------
print("Loading ...")
defm = pd.read_csv(DEF_FILE, dtype={"genome_id": str})
pro  = pd.read_csv(PRO_FILE, dtype=str)
bac  = pd.read_csv(BAC_FILE, dtype=str)

sys_cols = [c for c in defm.columns if c != "genome_id"]
meta = parse_def_columns(sys_cols)

# presence/absence (count > 0) indexed by host genome_id
counts = defm.set_index("genome_id")[sys_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
pres = (counts > 0).astype(int)

# collapse to the chosen granularity
if GRAIN == "type":
    grp = {}
    for c in sys_cols:
        grp.setdefault(meta[c]["type"], []).append(c)
    # present at TYPE level if ANY subtype present
    pres = pd.DataFrame({t: pres[cs].max(axis=1) for t, cs in grp.items()}, index=pres.index)
    sys_meta = pd.DataFrame({"system": list(grp.keys())})
    sys_meta["type"] = sys_meta["system"]
    sys_meta["subtypes"] = [",".join(sorted({meta[c]["subtype"] for c in grp[t]})) for t in grp]
    sys_meta["outcomes"] = [",".join(sorted({meta[c]["outcome"] for c in grp[t]})) for t in grp]
else:  # subtype grain
    lab = {c: f'{meta[c]["type"]}:{meta[c]["subtype"]}' for c in sys_cols}
    # Group by label to handle duplicate subtypes (combine with max/OR logic)
    grp = {}
    for c in sys_cols:
        grp.setdefault(lab[c], []).append(c)
    # present at SUBTYPE level if ANY occurrence present
    pres = pd.DataFrame({l: pres[cs].max(axis=1) for l, cs in grp.items()}, index=pres.index)
    # Build sys_meta with one row per unique subtype label
    sys_meta_rows = []
    for label in pres.columns:
        # Get first column in this group to extract metadata
        first_col = grp[label][0]
        sys_meta_rows.append({"system": label, "type": meta[first_col]["type"],
                              "subtype": meta[first_col]["subtype"], "outcome": meta[first_col]["outcome"]})
    sys_meta = pd.DataFrame(sys_meta_rows)

systems = list(pres.columns)
N_pop = pres.shape[0]                         # total S. suis genomes (background universe)
bg_carriers = pres.sum(axis=0)                # carriers per system across whole population
sys_meta = sys_meta.merge(
    (bg_carriers.rename("pop_carriers") / N_pop).rename("pop_freq").reset_index()
    .rename(columns={"index": "system"}), on="system", how="left")
sys_meta["pop_carriers"] = bg_carriers.reindex(sys_meta["system"]).values
sys_meta.to_csv(os.path.join(OUT_DIR, "defence_system_metadata.csv"), index=False)
print(f"  {N_pop} host genomes | {len(systems)} defence systems at '{GRAIN}' grain")

# ----------------------------------------------------------------------------
# 2. Build prophage -> host -> communities table
# ----------------------------------------------------------------------------
pro["prophage_id"] = pro["Genome_ID"]
pro["host"]        = pro["Genome_ID"].str.split("_", n=1).str[0]
pro["pro_comm"]    = "PRO_com_" + pro["Community"]

bac["host"]        = bac["Genome_ID"].str.replace(r"\.fna$", "", regex=True)
bac["bac_comm"]    = "BAC_com_" + bac["Community"]
host2bac = dict(zip(bac["host"], bac["bac_comm"]))

pro["bac_comm"] = pro["host"].map(host2bac)
long = pro[["prophage_id", "pro_comm", "host", "bac_comm"]].copy()
long.to_csv(os.path.join(OUT_DIR, "prophage_host_community_long.csv"), index=False)

# ----------------------------------------------------------------------------
# 3. Per prophage-community summary (spread across host communities)
# ----------------------------------------------------------------------------
summ = (long.groupby("pro_comm")
            .agg(n_prophages=("prophage_id", "size"),
                 n_distinct_hosts=("host", "nunique"),
                 n_bac_communities=("bac_comm", "nunique"))
            .reset_index())
bac_lists = (long.groupby("pro_comm")["bac_comm"]
                 .agg(lambda s: ",".join(sorted(s.dropna().unique())))
                 .rename("bac_communities"))
summ = summ.merge(bac_lists, on="pro_comm").sort_values(
    ["n_bac_communities", "n_distinct_hosts"], ascending=False)

# ----------------------------------------------------------------------------
# 4. Defence-system carriage among each community's distinct hosts
# ----------------------------------------------------------------------------
pairs = long[["pro_comm", "host"]].drop_duplicates()
merged = pairs.merge(pres, left_on="host", right_index=True, how="left")
carry  = merged.groupby("pro_comm")[systems].sum()                 # # hosts carrying each system
nhost  = pairs.groupby("pro_comm")["host"].nunique().rename("n_hosts")
frac   = carry.div(nhost, axis=0)                                  # fraction of hosts carrying

car_long = (carry.stack().rename("n_hosts_carrying").reset_index()
            .rename(columns={"level_1": "system"}))
car_long = car_long.merge(nhost.reset_index(), on="pro_comm")
car_long["frac_hosts_carrying"] = car_long["n_hosts_carrying"] / car_long["n_hosts"]
car_long["flag"] = np.select(
    [car_long["frac_hosts_carrying"] == 1.0,
     car_long["n_hosts_carrying"] == 0],
    ["always", "never"], default="variable")
car_long = car_long.merge(sys_meta[["system", "type", "pop_freq"]], on="system", how="left")
car_long.to_csv(os.path.join(OUT_DIR, "community_defence_carriage.csv"), index=False)

# ----------------------------------------------------------------------------
# 5. Enrichment / depletion test (binomial exact test)
#    For each (community, system) pair:
#      n = hosts in community
#      k = hosts carrying system
#      p = population frequency
#    Test: P(observing ≥ k carriers | n trials, probability p)
# ----------------------------------------------------------------------------
test_comms = nhost[nhost >= MIN_HOSTS].index
print(f"  testing {len(test_comms)} communities (>= {MIN_HOSTS} hosts) x {len(systems)} systems ...")
rows = []
for comm in test_comms:
    n = int(nhost[comm])
    for s in systems:
        k = int(carry.loc[comm, s])
        p = bg_carriers[s] / N_pop
        
        if k == 0 and bg_carriers[s] == 0:     # system absent everywhere -> uninformative
            continue
        
        # Binomial test: two-sided test for enrichment/depletion
        # p_value = P(X >= k | n, p) + P(X <= k | n, p) [two-sided]
        # For efficiency, use sf (survival function) and pmf (probability mass function)
        expected = n * p
        
        if k > expected:
            # Enriched: P(X >= k | n, p)
            p_value = binom.sf(k - 1, n, p)
            direction = "enriched"
        else:
            # Depleted: P(X <= k | n, p)
            p_value = binom.cdf(k, n, p)
            direction = "depleted"
        
        # Make it two-sided by doubling p-value (conservative)
        p_value = min(2 * p_value, 1.0)
        
        rows.append((comm, s, n, k, k / n, bg_carriers[s], p, p_value, direction))

enr = pd.DataFrame(rows, columns=[
    "pro_comm", "system", "n_hosts", "n_hosts_carrying", "frac_hosts_carrying",
    "pop_carriers", "pop_freq", "p_value", "direction"])
enr["q_value"] = benjamini_hochberg(enr["p_value"].values)
enr["significant"] = enr["q_value"] < FDR_ALPHA
enr = enr.merge(sys_meta[["system", "type"]], on="system", how="left")
enr = enr.sort_values(["significant", "q_value", "pro_comm"], ascending=[False, True, True])
enr.to_csv(os.path.join(OUT_DIR, "community_defence_enrichment.csv"), index=False)

# attach a quick "always / never" digest to the community summary
def digest(comm, want):
    sub = car_long[(car_long.pro_comm == comm) & (car_long.flag == want)]
    # for 'never', only list systems that exist in the population (informative absences)
    if want == "never":
        sub = sub[sub.pop_freq > 0]
    return ",".join(sorted(sub["system"]))
summ["always_systems"] = summ["pro_comm"].map(lambda c: digest(c, "always"))
n_sig = (enr[enr.significant].groupby("pro_comm").size().rename("n_sig_systems"))
summ = summ.merge(n_sig, on="pro_comm", how="left")
summ["n_sig_systems"] = summ["n_sig_systems"].fillna(0).astype(int)
summ.to_csv(os.path.join(OUT_DIR, "prophage_community_summary.csv"), index=False)

# ----------------------------------------------------------------------------
# 6. Console summary
# ----------------------------------------------------------------------------
print("\n--- Prophage communities spanning the most bacterial communities ---")
print(summ.head(12)[["pro_comm", "n_prophages", "n_distinct_hosts",
                     "n_bac_communities", "n_sig_systems"]].to_string(index=False))
sig = enr[enr.significant]
print(f"\nSignificant (q<{FDR_ALPHA}) community-system associations: {len(sig)} "
      f"({(sig.direction=='enriched').sum()} enriched, {(sig.direction=='depleted').sum()} depleted)")
print("\nTop enrichments (multi-host-community prophage communities):")
multi = set(summ[summ.n_bac_communities > 1]["pro_comm"])
top = sig[sig.pro_comm.isin(multi)].sort_values("q_value").head(15)
print(top[["pro_comm", "system", "n_hosts", "frac_hosts_carrying",
           "pop_freq", "direction", "q_value"]].to_string(index=False))

print("\nOutputs written to", os.path.abspath(OUT_DIR))
for f in ["prophage_host_community_long.csv", "prophage_community_summary.csv",
          "community_defence_carriage.csv", "community_defence_enrichment.csv",
          "defence_system_metadata.csv"]:
    print("  -", f)
