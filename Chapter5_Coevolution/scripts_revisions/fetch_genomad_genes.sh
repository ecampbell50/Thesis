#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Collect geNomad gene-level annotation for the genomes needed to answer
# supervisor comments 129, 143 and 191 on the coevolution chapter.
#
# RUN THIS ON THE HPC, from inside the 04_genomad directory
# (the one containing <genome>_genomad/ subdirectories).
#
#   1. copy hpc_genome_ids.txt into that directory
#   2. bash fetch_genomad_genes.sh
#   3. scp the resulting genomad_genes_for_ch5.tar.gz back to your laptop
#
# Collects, for each of the 78 genomes listed:
#   <id>_find_proviruses/<id>_provirus_genes.tsv     per-gene annotation of proviruses
#   <id>_find_proviruses/<id>_provirus_summary.tsv   provirus coordinates + scores
#   <id>_summary/<id>_virus_genes.tsv                per-gene annotation, summary set
#   <id>_summary/<id>_virus_summary.tsv              virus calls + taxonomy
# and, for the two medoids used in the alignment figure, the protein FASTA.
# ---------------------------------------------------------------------------
set -uo pipefail

IDS="${1:-hpc_genome_ids.txt}"
OUT="genomad_genes_for_ch5"
MEDOIDS="1307.2545 1307.3942"     # ICE_community_99 and PRO_com_74 representatives

[[ -f "$IDS" ]] || { echo "ERROR: cannot find $IDS in $(pwd)"; exit 1; }

mkdir -p "$OUT"
found=0; missing=0
: > "$OUT/_MISSING.txt"

while read -r g; do
  [[ -z "$g" ]] && continue
  base="${g}_genomad"
  if [[ ! -d "$base" ]]; then
    echo "$g  (no ${base}/ directory)" >> "$OUT/_MISSING.txt"; ((missing++)); continue
  fi
  got=0
  for rel in "${g}_find_proviruses/${g}_provirus_genes.tsv" \
             "${g}_find_proviruses/${g}_provirus_summary.tsv" \
             "${g}_summary/${g}_virus_genes.tsv" \
             "${g}_summary/${g}_virus_summary.tsv"; do
    if [[ -f "$base/$rel" ]]; then cp "$base/$rel" "$OUT/"; got=1; fi
  done
  if [[ $got -eq 1 ]]; then ((found++)); else
    echo "$g  (directory present, no expected files)" >> "$OUT/_MISSING.txt"; ((missing++))
  fi
done < "$IDS"

# protein sequences for the two figure medoids only (keeps the archive small)
for g in $MEDOIDS; do
  for rel in "${g}_genomad/${g}_summary/${g}_virus_proteins.faa" \
             "${g}_genomad/${g}_find_proviruses/${g}_provirus_proteins.faa"; do
    [[ -f "$rel" ]] && cp "$rel" "$OUT/"
  done
done

printf '%s\n' "genomes with files collected : $found" \
              "genomes missing              : $missing" \
              "files collected              : $(ls -1 "$OUT" | grep -vc '^_MISSING' || true)" \
  | tee "$OUT/_SUMMARY.txt"
[[ -s "$OUT/_MISSING.txt" ]] && echo "(see $OUT/_MISSING.txt)" || rm -f "$OUT/_MISSING.txt"

tar -czf "${OUT}.tar.gz" "$OUT"
echo
echo "Done -> ${OUT}.tar.gz  ($(du -h "${OUT}.tar.gz" | cut -f1))"
echo "Copy it back with, from your laptop:"
echo "  scp <user>@login2:\$(pwd)/${OUT}.tar.gz ~/Desktop/Thesis/1_CoevolutionofSsuisProphages/"
