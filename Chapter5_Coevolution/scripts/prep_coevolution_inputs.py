#!/usr/bin/env python3
# Host x prophage community incidence matrix (data/host_by_procomm.csv).
# Prophage communities with >= 5 distinct hosts only, as in plot_bipartite_network.R.

import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # Chapter5_Coevolution/
MIN_COMM_HOSTS = 5

long = pd.read_csv(os.path.join(ROOT, "data/prophage_host_community_long.csv"),
                   dtype={"host": str, "pro_comm": str, "bac_comm": str})

comm_hostcount = long.groupby("pro_comm")["host"].nunique()
real_comms = comm_hostcount[comm_hostcount >= MIN_COMM_HOSTS].index

sub = long[long["pro_comm"].isin(real_comms)]
host_by_procomm = (pd.crosstab(sub["host"], sub["pro_comm"]) > 0).astype(int)
host_by_procomm.to_csv(os.path.join(ROOT, "data/host_by_procomm.csv"))
