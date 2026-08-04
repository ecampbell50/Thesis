#!/usr/bin/env bash
# Is the PRO_com_74 representative contained within the ICE_community_99 representative?
# Follow-up to VIRIDIC, which scored these two at 50.1% intergenomic similarity.
# Chapter 5. Written 2026-07-30.
set -euo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)/Top15MedoidProphages"
ICE="$DIR/1307.2545_PKOW01000001_provirus-114-105135.fna"        # ICE_community_99 medoid, 105,022 bp
PHAGE="$DIR/1307.3942_JAIMDQ010000003_provirus-66207-113519.fna" # PRO_com_74 medoid,      47,313 bp
OUT="$(cd "$(dirname "$0")" && pwd)/pc74_vs_ice99"

makeblastdb -in "$ICE" -dbtype nucl -out "${OUT}_db" >/dev/null
blastn -query "$PHAGE" -db "${OUT}_db" -evalue 1e-5 \
       -outfmt "6 qstart qend sstart send pident length" \
       > "${OUT}_hits.tsv"

# merge overlapping query intervals to get non-redundant covered length
python3 - "$PHAGE" "${OUT}_hits.tsv" <<'PY'
import sys
qlen = sum(len(l.strip()) for l in open(sys.argv[1]) if not l.startswith('>'))
iv = []
for line in open(sys.argv[2]):
    f = line.split()
    a, b = int(f[0]), int(f[1]); iv.append((min(a, b), max(a, b)))
iv.sort(); merged = []
for s, e in iv:
    if merged and s <= merged[-1][1] + 1:
        merged[-1][1] = max(merged[-1][1], e)
    else:
        merged.append([s, e])
cov = sum(e - s + 1 for s, e in merged)
print(f'query length      : {qlen:,} bp')
print(f'aligned (merged)  : {cov:,} bp = {100*cov/qlen:.1f}%')
print(f'blocks            : {len(merged)}')
PY
