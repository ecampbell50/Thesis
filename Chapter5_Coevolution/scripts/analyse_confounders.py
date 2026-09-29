#!/usr/bin/env python3
# =============================================================================
# Q1. Does prophage burden predict host DEFENCE burden once confounders are
#     controlled (genome size, serotype/clade, assembly fragmentation)?
#     -> negative-binomial regression + partial Spearman sensitivity.
# Q2. Is prophage SIZE (provirus length, bp) associated with:
#     A) prophage community,  B) host defence burden,  C) community promiscuity?
#
# Inputs (all local):
#   data/genome_defence_by_location.csv  (per-genome defence by location + burden)
#   SupplementaryTable1-BVBRC_29Nov23.csv                       (Size, CDS, Contigs, N50, GC) - strip 'Ssuis_'
#   ConsensusSerotypes_23Jul24.csv                              (Serotype_Consensus) - strip '.gff'
#   data/prophage_host_community_long.csv   (prophage_id, pro_comm, host, bac_comm)
#   data/genome_community_mapping_44_PRO.csv                 (provirus coords in the ID -> length)
# =============================================================================

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, kruskal, rankdata
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")

OUT_DIR = "tables"
os.makedirs(OUT_DIR, exist_ok=True)


def partial_spearman(df, x, y, covars):
    """Spearman partial correlation of x,y controlling for covars (rank-residual method)."""
    d = df[[x, y] + covars].dropna()
    R = {c: rankdata(d[c]) for c in [x, y] + covars}
    X = np.column_stack([np.ones(len(d))] + [R[c] for c in covars])
    def resid(v):
        beta, *_ = np.linalg.lstsq(X, v, rcond=None)
        return v - X @ beta
    rx, ry = resid(R[x]), resid(R[y])
    r = np.corrcoef(rx, ry)[0, 1]
    n = len(d); k = len(covars)
    from scipy.stats import t as tdist
    tval = r * np.sqrt((n - k - 2) / (1 - r**2))
    p = 2 * tdist.sf(abs(tval), n - k - 2)
    return r, p, n


# ----------------------------------------------------------------------------
# Load + build master per-genome table
# ----------------------------------------------------------------------------
print("Loading ...")
g = pd.read_csv("data/genome_defence_by_location.csv", dtype={"genome": str})
bv = pd.read_csv("data/SupplementaryTable1-BVBRC_29Nov23.csv", dtype=str)
bv["genome"] = bv["Genome ID"].str.replace(r"^Ssuis_", "", regex=True)
sero = pd.read_csv("../Chapter3_DefenceProfile/data/ConsensusSerotypes_23Jul24.csv", dtype=str)
sero["genome"] = sero["Strain"].str.replace(r"\.gff$", "", regex=True)

m = (g.merge(bv[["genome", "Size", "CDS", "Contigs", "Contig N50", "GC Content"]], on="genome", how="left")
       .merge(sero[["genome", "Serotype_Consensus"]], on="genome", how="left"))
for c in ["Size", "CDS", "Contigs", "Contig N50", "GC Content"]:
    m[c] = pd.to_numeric(m[c], errors="coerce")
m["log_size"]    = np.log10(m["Size"])
m["log_contigs"] = np.log10(m["Contigs"])
m["gc"]          = m["GC Content"]
# pool rare serotypes / untypable into 'other'; reference = serotype 2 (largest)
keep = m.loc[m.Serotype_Consensus != "-", "Serotype_Consensus"].value_counts()
keep = set(keep[keep >= 20].index)
m["sero"] = np.where(m.Serotype_Consensus.isin(keep), m.Serotype_Consensus, "other")
m = m.dropna(subset=["log_size", "gc", "n_def_host"]).copy()
print(f"  {len(m)} genomes modelled | {m.sero.nunique()} serotype levels (ref=2)")

# ============================================================================
# Q1. Multivariable model: host-core defence ~ burden + confounders
# ============================================================================
print("\n" + "=" * 70 + "\nQ1. BURDEN -> HOST-CORE DEFENCE, adjusted\n" + "=" * 70)

print("\nStaged Spearman (showing where the raw signal goes):")
for label, ycol, cov in [
    ("raw: burden vs TOTAL defence",            "n_def_genes", []),
    ("burden vs HOST-CORE defence",             "n_def_host",  []),
    ("burden vs host-core | genome size",       "n_def_host",  ["log_size"]),
    ("burden vs host-core | size + fragmentation", "n_def_host", ["log_size", "log_contigs"])]:
    if cov:
        r, p, n = partial_spearman(m, "n_prophages", ycol, cov)
        print(f"  {label:46s} partial rho={r:+.3f}  p={p:.1e}")
    else:
        r, p = spearmanr(m["n_prophages"], m[ycol])
        print(f"  {label:46s}         rho={r:+.3f}  p={p:.1e}")

# Negative-binomial regression (estimates dispersion). exp(coef) = IRR.
def run_nb(formula, data, tag):
    mod = smf.negativebinomial(formula, data=data).fit(disp=0, maxiter=200)
    irr = np.exp(mod.params["n_prophages"])
    ci = np.exp(mod.conf_int().loc["n_prophages"])
    print(f"\n  [{tag}]  burden IRR = {irr:.4f}  (95% CI {ci[0]:.4f}-{ci[1]:.4f})  "
          f"p = {mod.pvalues['n_prophages']:.2e}")
    return mod

print("\nNegative-binomial models (outcome = host-core defence gene count):")
m_uni  = run_nb("n_def_host ~ n_prophages", m, "unadjusted")
m_size = run_nb("n_def_host ~ n_prophages + log_size", m, "+ genome size")
m_full = run_nb("n_def_host ~ n_prophages + log_size + log_contigs + gc + C(sero, Treatment('2'))",
                m, "+ size + fragmentation + GC + serotype")
# sensitivity: swap genome size for CDS
m_cds  = run_nb("n_def_host ~ n_prophages + np.log10(CDS) + log_contigs + gc + C(sero, Treatment('2'))",
                m, "sensitivity: CDS instead of Size")

# tidy coefficient table from the full model
coef = pd.DataFrame({"term": m_full.params.index, "coef": m_full.params.values,
                     "IRR": np.exp(m_full.params.values), "p_value": m_full.pvalues.values})
coef = coef.merge(np.exp(m_full.conf_int()).rename(columns={0: "IRR_lo", 1: "IRR_hi"}),
                  left_on="term", right_index=True)
coef.to_csv(os.path.join(OUT_DIR, "Q1_nb_full_model_coefficients.csv"), index=False)
m.to_csv(os.path.join(OUT_DIR, "genome_master_table.csv"), index=False)

# ============================================================================
# Q2. Prophage SIZE associations
# ============================================================================
print("\n" + "=" * 70 + "\nQ2. PROPHAGE SIZE associations\n" + "=" * 70)
pro  = pd.read_csv("data/genome_community_mapping_44_PRO.csv", dtype=str)
long = pd.read_csv("data/prophage_host_community_long.csv", dtype=str)
ext = pro["Genome_ID"].str.extract(r"_provirus-(?P<s>\d+)-(?P<e>\d+)\.fna$").astype(float)
pro["length"] = ext["e"] - ext["s"] + 1
pro = pro.rename(columns={"Genome_ID": "prophage_id"})[["prophage_id", "length"]]
P = long.merge(pro, on="prophage_id", how="left").dropna(subset=["length"])
# community promiscuity = # distinct bac communities a prophage community connects to
promisc = P.dropna(subset=["bac_comm"]).groupby("pro_comm")["bac_comm"].nunique().rename("n_bac_comm")
P = P.merge(promisc, on="pro_comm", how="left")
P["promisc_class"] = np.where(P["n_bac_comm"] > 1, "promiscuous", "single")

# --- A) length vs prophage community ---
comm_n = P.groupby("pro_comm")["length"].size()
big = comm_n[comm_n >= 10].index
samples = [P.loc[P.pro_comm == c, "length"].values for c in big]
H, p = kruskal(*samples)
# variance explained by community (eta^2 on ranks ~ Kruskal effect size)
eta2 = (H - len(big) + 1) / (len(P[P.pro_comm.isin(big)]) - len(big))
print(f"\nA) length ~ prophage community: Kruskal-Wallis H={H:.0f}, p={p:.1e} "
      f"across {len(big)} communities (n>=10); rank eta^2={eta2:.2f}")

# --- B) length vs host defence burden (per genome) ---
perhost = (P.groupby("host").agg(total_prophage_len=("length", "sum"),
                                 mean_prophage_len=("length", "mean"),
                                 n_prophages=("prophage_id", "size")).reset_index())
perhost = perhost.merge(m[["genome", "n_def_host", "n_def_genes"]],
                        left_on="host", right_on="genome", how="inner")
print("\nB) prophage size vs host defence (per genome, prophage-carrying hosts):")
for xcol in ["total_prophage_len", "mean_prophage_len"]:
    for ycol in ["n_def_host", "n_def_genes"]:
        r, pp = spearmanr(perhost[xcol], perhost[ycol])
        print(f"   {xcol:20s} vs {ycol:12s} rho={r:+.3f}  p={pp:.1e}")
# does mean prophage size add to defence beyond count? partial controlling for count + size
mh = perhost.merge(m[["genome", "log_size"]], on="genome", how="left")
r, pp, n = partial_spearman(mh, "mean_prophage_len", "n_def_host", ["n_prophages", "log_size"])
print(f"   mean_prophage_len vs n_def_host | (count + genome size): partial rho={r:+.3f} p={pp:.1e}")

# --- C) length vs community promiscuity ---
u, pmw = mannwhitneyu(P.loc[P.promisc_class == "promiscuous", "length"],
                      P.loc[P.promisc_class == "single", "length"], alternative="two-sided")
med_p = P.loc[P.promisc_class == "promiscuous", "length"].median()
med_s = P.loc[P.promisc_class == "single", "length"].median()
print(f"\nC) prophage length vs promiscuity (per prophage): "
      f"promiscuous median={med_p:.0f} vs single median={med_s:.0f}  MWU p={pmw:.1e}")
comm_lvl = P.groupby("pro_comm").agg(mean_len=("length", "mean"), n_bac_comm=("n_bac_comm", "first"))
rc, pc = spearmanr(comm_lvl["mean_len"], comm_lvl["n_bac_comm"])
print(f"   community-level: mean length vs n_bac_communities  rho={rc:+.3f}  p={pc:.1e}")
P.to_csv(os.path.join(OUT_DIR, "prophage_length_table.csv"), index=False)

# ============================================================================
# Figures
# ============================================================================
# Q1: staged attenuation of the burden effect (Spearman)
stages, vals = [], []
stages.append("raw\n(total defence)"); vals.append(spearmanr(m.n_prophages, m.n_def_genes)[0])
stages.append("host-core only");      vals.append(spearmanr(m.n_prophages, m.n_def_host)[0])
vals.append(partial_spearman(m, "n_prophages", "n_def_host", ["log_size"])[0]); stages.append("| genome size")
vals.append(partial_spearman(m, "n_prophages", "n_def_host", ["log_size", "log_contigs"])[0]); stages.append("| size+frag")
fig, ax = plt.subplots(figsize=(7, 5))
ax.bar(stages, vals, color=["#6baed6", "#4292c6", "#2171b5", "#08519c"])
ax.axhline(0, color="k", lw=.8)
ax.set_ylabel("Spearman rho (burden vs defence)")
ax.set_title("Burden->defence association attenuates as confounders are removed")
for i, v in enumerate(vals):
    ax.text(i, v + 0.005, f"{v:.2f}", ha="center", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, "fig_Q1_attenuation.png"), dpi=200); plt.close(fig)

# Q2C: prophage length by promiscuity
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
sns.boxplot(data=P, x="promisc_class", y="length", order=["single", "promiscuous"],
            hue="promisc_class", palette={"single": "#bdbdbd", "promiscuous": "#08519c"},
            legend=False, showfliers=False, ax=ax[0])
ax[0].set_yscale("log"); ax[0].set_ylabel("provirus length (bp, log)")
ax[0].set_xlabel("prophage community type"); ax[0].set_title("C) prophage size vs promiscuity")
# Q2A: length by the largest communities
top_comms = comm_n.sort_values(ascending=False).head(15).index
sub = P[P.pro_comm.isin(top_comms)]
order_c = sub.groupby("pro_comm")["length"].median().sort_values().index
sns.boxplot(data=sub, x="pro_comm", y="length", order=order_c, showfliers=False,
            color="#41ab5d", ax=ax[1])
ax[1].set_yscale("log"); ax[1].set_ylabel("provirus length (bp, log)")
ax[1].set_xlabel("prophage community"); ax[1].set_title("A) size differs by community")
ax[1].tick_params(axis="x", rotation=90)
fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, "fig_Q2_prophage_size.png"), dpi=200); plt.close(fig)

print("\nOutputs ->", os.path.abspath(OUT_DIR))
for f in ["genome_master_table.csv", "Q1_nb_full_model_coefficients.csv",
          "prophage_length_table.csv", "fig_Q1_attenuation.png", "fig_Q2_prophage_size.png"]:
    print("  -", f)
