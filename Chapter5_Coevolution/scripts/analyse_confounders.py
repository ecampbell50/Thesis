#!/usr/bin/env python3
# Per-genome and per-prophage tables for BipartiteNetwork_Analysis.ipynb (Fig S5.6, S5.7).
# Writes data/confounders/genome_master_table.csv and prophage_length_table.csv.

import os
import numpy as np
import pandas as pd

OUT_DIR = "data/confounders"
os.makedirs(OUT_DIR, exist_ok=True)

# Per-genome table: defence counts, assembly stats, serotype
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
# serotypes with < 20 genomes (and untypable) pooled as 'other'
keep = m.loc[m.Serotype_Consensus != "-", "Serotype_Consensus"].value_counts()
keep = set(keep[keep >= 20].index)
m["sero"] = np.where(m.Serotype_Consensus.isin(keep), m.Serotype_Consensus, "other")
m = m.dropna(subset=["log_size", "gc", "n_def_host"]).copy()
m.to_csv(os.path.join(OUT_DIR, "genome_master_table.csv"), index=False)

# Per-prophage table: length (from the geNomad coordinates in the ID) and community promiscuity
pro  = pd.read_csv("data/genome_community_mapping_44_PRO.csv", dtype=str)
long = pd.read_csv("data/prophage_host_community_long.csv", dtype=str)
ext = pro["Genome_ID"].str.extract(r"_provirus-(?P<s>\d+)-(?P<e>\d+)\.fna$").astype(float)
pro["length"] = ext["e"] - ext["s"] + 1
pro = pro.rename(columns={"Genome_ID": "prophage_id"})[["prophage_id", "length"]]
P = long.merge(pro, on="prophage_id", how="left").dropna(subset=["length"])
promisc = P.dropna(subset=["bac_comm"]).groupby("pro_comm")["bac_comm"].nunique().rename("n_bac_comm")
P = P.merge(promisc, on="pro_comm", how="left")
P["promisc_class"] = np.where(P["n_bac_comm"] > 1, "promiscuous", "single")
P.to_csv(os.path.join(OUT_DIR, "prophage_length_table.csv"), index=False)
