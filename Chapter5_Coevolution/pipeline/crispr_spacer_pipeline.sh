#!/bin/bash
# ============================================================================
# crispr_spacer_pipeline.sh  —  Move 3 (HPC): CRISPR spacer -> prophage matching
#
# A host CRISPR spacer that matches a prophage in the dataset is a DIRECT,
# recorded coevolutionary event: "this lineage was immunised against this phage."
# This overlays an immunity record on the bipartite network.
#
# Run on Kelvin2 (adjust the SBATCH header / paths if your layout differs).
#   sbatch crispr_spacer_pipeline.sh
#
# Expects (from your existing layout, see Methodology.md):
#   $BASE/01_genomes/*.fna        2,119 host assemblies
#   $BASE/05_proviruses/*.fna     3,216 extracted proviruses
# Produces:
#   $BASE/09_crispr/spacers.fasta                 all spacers (header: >GENOME|array|spacerN)
#   $BASE/09_crispr/spacer_vs_provirus.tsv        filtered blastn hits  <-- pull this locally
# ============================================================================
#SBATCH --job-name=crispr_spacers
#SBATCH --partition=k2-bioinf,k2-lowpri
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16G
#SBATCH --output=crispr_spacers.out
#SBATCH --error=crispr_spacers.err
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=END,FAIL

set -euo pipefail
BASE=/mnt/scratch2/users/40204129/Chapter3_Coevolution
GENOMES=$BASE/01_genomes
PROVIR=$BASE/05_proviruses
OUT=$BASE/09_crispr
mkdir -p "$OUT"/spacers_raw

# --- env: minced is already in the `prokka` env; blast (makeblastdb/blastn) ships
#     with prokka too. If seqkit is missing, the awk relabel in step 3 covers it.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate prokka

# ---------------------------------------------------------------------------
# 1. Detect CRISPR arrays per genome (minced), then build the tagged spacer FASTA.
#    CONFIRMED behaviour (minced 0.4.2): the spacer FASTA is written NEXT TO the
#    .crisprs output (2nd positional) as "<genome-basename>_spacers.fa". Many will
#    be empty (genomes with no array). Headers carry the CONTIG name, not the genome
#    id, so step 1b prepends "${gid}|" for the downstream join. Glob must be *.fna.
#    minced defaults (minSL 26 / maxSL 50 / minNR 3) are standard for S. suis.
# ---------------------------------------------------------------------------
echo "[1a] minced on $(ls "$GENOMES"/*.fna | wc -l) genomes ..."
mkdir -p "$OUT/arrays"
for f in "$GENOMES"/*.fna; do
    gid=$(basename "$f" .fna)
    minced -gffFull -spacers "$f" "$OUT/arrays/${gid}.crisprs" "$OUT/arrays/${gid}.gff" \
        >/dev/null 2>&1 || true       # spacers land at $OUT/arrays/${gid}_spacers.fa
done

echo "[1b] building genome-tagged all_spacers.fasta ..."
: > "$OUT/spacers.fasta"
for sp in "$OUT"/arrays/*_spacers.fa; do
    [[ -s "$sp" ]] || continue                       # skip genomes with no arrays
    gid=$(basename "$sp" _spacers.fa)
    sed "s/^>/>${gid}|/" "$sp" >> "$OUT/spacers.fasta"
done
echo "    spacers: $(grep -c '^>' "$OUT/spacers.fasta") from $(grep '^>' "$OUT/spacers.fasta" | cut -d'|' -f1 | sort -u | wc -l) genomes"

# ---------------------------------------------------------------------------
# 2. (optional but recommended) keep only spacers from genomes DefenseFinder
#    flagged CRISPR-Cas positive -> removes spurious minced arrays.
#    Provide a one-column file of CRISPR+ genome ids at $OUT/crispr_pos_genomes.txt
# ---------------------------------------------------------------------------
if [[ -s "$OUT/crispr_pos_genomes.txt" ]]; then
    echo "[2] filtering spacers to DefenseFinder CRISPR+ genomes ..."
    seqkit grep -r -n -f <(sed 's/^/^/' "$OUT/crispr_pos_genomes.txt") "$OUT/spacers.fasta" \
        > "$OUT/spacers.filtered.fasta" || cp "$OUT/spacers.fasta" "$OUT/spacers.filtered.fasta"
else
    cp "$OUT/spacers.fasta" "$OUT/spacers.filtered.fasta"
fi

# ---------------------------------------------------------------------------
# 3. Provirus BLAST database (all 3,216 proviruses concatenated)
# ---------------------------------------------------------------------------
echo "[3] building provirus blast db ..."
# Relabel every provirus sequence header to its FILENAME (with .fna) so that
# blast sseqid == prophage_id in prophage_host_community_long.csv -> trivial join.
: > "$OUT/all_proviruses.fna"
for f in "$PROVIR"/*.fna; do
    id=$(basename "$f")          # e.g. 1004951.3_CP002651_provirus-1291910-1341951.fna
    # rewrite every header line to the filename (awk; no seqkit dependency)
    awk -v id="$id" '/^>/{print ">"id; next}{print}' "$f" >> "$OUT/all_proviruses.fna"
done
makeblastdb -in "$OUT/all_proviruses.fna" -dbtype nucl -out "$OUT/provirus_db" >/dev/null

# ---------------------------------------------------------------------------
# 4. Match spacers -> proviruses.  Short query => blastn-short, dust off.
# ---------------------------------------------------------------------------
echo "[4] blastn spacers vs proviruses ..."
blastn -task blastn-short -query "$OUT/spacers.filtered.fasta" -db "$OUT/provirus_db" \
    -dust no -evalue 1 -num_threads "${SLURM_NTASKS:-8}" \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore' \
    > "$OUT/spacer_vs_provirus.raw.tsv"

# ---------------------------------------------------------------------------
# 5. Protospacer-grade filter: near-full-length, high identity, <=1 mismatch.
#    (cols: 1 qseqid 2 sseqid 3 pident 4 length 5 mismatch 6 gapopen 7 qlen ...)
# ---------------------------------------------------------------------------
echo "[5] filtering hits (pident>=95, mismatch<=1, coverage>=0.9) ..."
awk -F'\t' '($3>=95) && ($5<=1) && ($4/$7>=0.90)' "$OUT/spacer_vs_provirus.raw.tsv" \
    > "$OUT/spacer_vs_provirus.tsv"
echo "    kept $(wc -l < "$OUT/spacer_vs_provirus.tsv") protospacer-grade hits"

echo "DONE. Pull with:"
echo "  scp -i ~/.ssh/my-kelvin-key 40204129@kelvin2.qub.ac.uk:$OUT/spacer_vs_provirus.tsv \\"
echo "      \"$BASE_LOCAL/Results/coevolution/\"   # then run analyse_crispr_targeting.py"
