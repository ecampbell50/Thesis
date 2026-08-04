#!/usr/bin/env python3
"""
Patristic distance summaries used in Chapter 5:
  - clonality of bacterial community 2 (median / max within-clade distance)
  - species-wide median and tree diameter for scale
  - carriage of prophage community 4 by bacterial community

Requires treedist_all (Creevey, doi:10.5281/zenodo.1244019), see ../treedist/.
Run first:  ./treedist/treedist_all <tree> matrix     -> <tree>.dist
Written 2026-07-30.
"""
import csv, sys, statistics, itertools
from collections import Counter

DIST = sys.argv[1] if len(sys.argv) > 1 else 'SsuisPhylo_FastTree_11Sep24.treefile.dist'
R = '../Results/'

fh = open(DIST)
names = fh.readline().rstrip('\n').split('\t')[1:]
idx = {n: i for i, n in enumerate(names)}
D = {}
for line in fh:
    p = line.rstrip('\n').rstrip('\t').split('\t')
    if len(p) > 1:
        D[p[0]] = [float(x) for x in p[1:]]
d = lambda a, b: D[a][idx[b]]

com = {r['Genome_ID'].replace('.fna', ''): r['Community']
       for r in csv.DictReader(open(R + 'genome_community_mapping_49_BAC.csv'))}
p4 = {r['host'] for r in csv.DictReader(open(R + 'confounders/prophage_length_table.csv'))
      if r['pro_comm'] == 'PRO_com_4'}

bac2 = [g for g, c in com.items() if c == '2' and g in idx]
w = [d(a, b) for a, b in itertools.combinations(bac2, 2)]
print(f'bacterial community 2: n={len(bac2)}  median={statistics.median(w):.5f}  max={max(w):.5f}')

ss = [n for n in names if not n.startswith('Spneumoniae')]   # exclude outgroup
allp = [D[a][idx[b]] for i, a in enumerate(ss) for b in ss[i+1:]]
print(f'all S. suis pairs: n={len(allp):,}  median={statistics.median(allp):.5f}  '
      f'mean={statistics.mean(allp):.5f}  diameter={max(allp):.5f}')

sizes = Counter(com.values())
print('\nPRO_com_4 carriage by bacterial community:')
for c in sorted({com[h] for h in p4 if h in com}, key=lambda x: -sizes[x]):
    n = sum(1 for h in p4 if com.get(h) == c)
    print(f'  community {c:>4}: {n:>4} / {sizes[c]:<4} ({100*n/sizes[c]:.1f}%)')
