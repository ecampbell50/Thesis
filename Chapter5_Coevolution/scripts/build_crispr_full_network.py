#!/usr/bin/env python3
"""
build_crispr_full_network.py
Builds the FULL CRISPR tripartite network, made tractable two ways:
  * spacers  -> clustered by IDENTICAL sequence (one representative square each)
  * prophages-> collapsed to their PRO community (one triangle each)

Nodes: genome (circle) | prophage_community (triangle) | spacer_cluster (square)
Edges: genome-spacer_cluster  "spacer in genome"
       spacer_cluster-community "protospacer hit"
       genome-community         "prophage in genome"  (crimson if the genome ALSO
                                  has a spacer targeting that community = coexistence)

Run AFTER scp'ing all_spacers.fasta into data/.
Out: data/crispr_full_{nodes,edges}.csv  -> plot_crispr_full_network.R
"""
import os, sys
import pandas as pd

COE = "data"
FASTA = f"{COE}/all_spacers.fasta"
if not os.path.exists(FASTA):
    sys.exit(f"missing {FASTA}\nscp it from HPC: 09_minced/all_spacers.fasta")

# ---- 1. parse fasta, cluster identical spacer sequences -------------------
qid2seq = {}
name, seq = None, []
with open(FASTA) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None: qid2seq[name] = "".join(seq)
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name is not None: qid2seq[name] = "".join(seq)

sp = pd.DataFrame({"qseqid": list(qid2seq), "seq": list(qid2seq.values())})
sp["genome"] = sp.qseqid.str.split("|").str[0]

# spacer clusters: prefer cd-hit 90% FAMILIES if present, else exact-identity.
CLSTR = f"{COE}/spacers_cdhit90.clstr"
if os.path.exists(CLSTR):
    fam, cur = {}, None
    for line in open(CLSTR):
        if line.startswith(">Cluster"):
            cur = "FAM" + line.split()[1]
        elif ">" in line:
            nm = line.split(">", 1)[1].split("...")[0].strip()
            fam[nm] = cur
    sp["spacer_cluster"] = sp.qseqid.map(fam)
    method = "cd-hit 90% families"
else:
    sp["spacer_cluster"] = sp.groupby("seq").ngroup().map(lambda i: f"SPC{i:04d}")
    method = "exact-identity"
n_unique = sp.spacer_cluster.nunique()
print(f"spacer clustering: {method}")

# ---- 2. hits, mapped to PRO communities + spacer clusters ------------------
cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
        "qlen","qstart","qend","sstart","send","evalue","bitscore"]
h = pd.read_csv(f"{COE}/spacer_vs_provirus.tsv", sep="\t", names=cols, dtype=str)
long = pd.read_csv("data/prophage_host_community_long.csv", dtype=str)
prov2comm = dict(zip(long.prophage_id, long.pro_comm))
qid2clu = dict(zip(sp.qseqid, sp.spacer_cluster))
h["pro_comm"] = h.sseqid.map(prov2comm)
h["spacer_cluster"] = h.qseqid.map(qid2clu)
h = h.dropna(subset=["pro_comm", "spacer_cluster"])

# ---- 3. node + edge sets ---------------------------------------------------
hit_clusters = set(h.spacer_cluster)
hit_comms    = set(h.pro_comm)
genomes      = set(h.qseqid.str.split("|").str[0])             # genomes with a hitting spacer

sp_hit     = sp[sp.spacer_cluster.isin(hit_clusters) & sp.genome.isin(genomes)]
gen_spacer = sp_hit[["genome","spacer_cluster"]].drop_duplicates()
hit_edges  = h[["spacer_cluster","pro_comm"]].drop_duplicates()
res        = long[long.host.isin(genomes) & long.pro_comm.isin(hit_comms)][["host","pro_comm"]].drop_duplicates()
res.columns = ["genome","pro_comm"]

# coexistence: genome carries a community AND has a spacer cluster hitting it
gen_hits = gen_spacer.merge(hit_edges, on="spacer_cluster")[["genome","pro_comm"]].drop_duplicates()
gen_hits_set = set(map(tuple, gen_hits.values))
res["coexistence"] = [(g,c) in gen_hits_set for g,c in zip(res.genome, res.pro_comm)]

# ---- 4. write nodes + edges ------------------------------------------------
clu_ngen = sp_hit.groupby("spacer_cluster").genome.nunique()
nodes = pd.concat([
    pd.DataFrame({"name": sorted(genomes),      "type": "genome"}),
    pd.DataFrame({"name": sorted(hit_comms),    "type": "prophage_community"}),
    pd.DataFrame({"name": sorted(hit_clusters), "type": "spacer_cluster"}),
], ignore_index=True)
nodes["n_genomes"] = nodes.apply(
    lambda r: int(clu_ngen.get(r["name"],1)) if r.type=="spacer_cluster" else 1, axis=1)

e1 = gen_spacer.rename(columns={"genome":"from","spacer_cluster":"to"}); e1["etype"]="spacer in genome";  e1["coexistence"]=False
e2 = hit_edges.rename(columns={"spacer_cluster":"from","pro_comm":"to"}); e2["etype"]="protospacer hit";   e2["coexistence"]=False
e3 = res.rename(columns={"genome":"from","pro_comm":"to"});               e3["etype"]="prophage in genome"; e3["coexistence"]=res.coexistence.values
E = pd.concat([e1,e2,e3], ignore_index=True)
E.loc[(E.etype=="prophage in genome") & (E.coexistence), "etype"] = "coexistence (carries a phage it targets)"

nodes.to_csv(f"{COE}/crispr_full_nodes.csv", index=False)
E.to_csv(f"{COE}/crispr_full_edges.csv", index=False)

print(f"spacers: {len(sp)} instances -> {n_unique} unique sequences ({len(hit_clusters)} with hits)")
print(f"nodes : {len(genomes)} genomes | {len(hit_comms)} prophage communities | {len(hit_clusters)} spacer clusters "
      f"= {len(nodes)} total")
print(f"edges : {len(e1)} spacer-in-genome | {len(e2)} hits | {len(e3)} resident "
      f"({int(res.coexistence.sum())} coexistence) = {len(E)} total")
print("wrote crispr_full_nodes.csv + crispr_full_edges.csv")
