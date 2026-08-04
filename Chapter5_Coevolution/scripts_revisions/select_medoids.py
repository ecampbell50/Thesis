#!/usr/bin/env python3
"""
Select a medoid representative for each of the 15 largest prophage communities.

The medoid is the community member with the highest MEAN within-community
Sourmash Jaccard similarity to all other members, i.e. the most typical
sequence of its community. Used as VIRIDIC input (Chapter 5, Section 5.2.x).

Inputs  (relative to ../Results/):
  confounders/prophage_length_table.csv                  prophage_id -> pro_comm
  Ssuis_3216proviruses_k15s500_44_intra_community_PRO.csv  Source,Target,Value (Jaccard)
Output:
  ../VIRIDIC_top15_representatives.txt
Written 2026-07-30.
"""
import csv, os
from collections import defaultdict

R = os.path.join(os.path.dirname(__file__), '..', 'Results')

members = defaultdict(list)
with open(os.path.join(R, 'confounders/prophage_length_table.csv')) as fh:
    for r in csv.DictReader(fh):
        members[r['pro_comm']].append(r['prophage_id'])

top = sorted(members.items(), key=lambda kv: -len(kv[1]))[:15]
topset = {c for c, _ in top}
memberof = {x: c for c, m in members.items() for x in m}

tot, cnt = defaultdict(float), defaultdict(int)
with open(os.path.join(R, 'Ssuis_3216proviruses_k15s500_44_intra_community_PRO.csv')) as fh:
    for r in csv.DictReader(fh):
        a, b, v = r['Source'], r['Target'], float(r['Value'])
        if a == b:
            continue
        ca, cb = memberof.get(a), memberof.get(b)
        if ca is None or ca != cb or ca not in topset:
            continue
        tot[a] += v; cnt[a] += 1
        tot[b] += v; cnt[b] += 1

for comm, mem in top:
    scored = sorted(((tot[x] / cnt[x], x) for x in mem if cnt.get(x)), reverse=True)
    if scored:
        print(f'{comm}\t{len(mem)}\t{scored[0][0]:.4f}\t{scored[0][1]}')
