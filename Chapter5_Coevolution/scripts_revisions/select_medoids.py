#!/usr/bin/env python3
# Medoid of each of the 15 largest prophage communities, for VIRIDIC (Fig 5.6).

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
