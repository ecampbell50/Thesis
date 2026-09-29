#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Figure 3.1 and Figure S3.1
# MMseqs2 clustering (95% id, 80% cov) + all-vs-all search that builds the network edges.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/DefGenes_Network/mmseqs2_80clustcov_16APR26.sh
# ======================================================================
#SBATCH --job-name=mmseqs2_80clustcov
#SBATCH --error=mmseqs2_80clustcov_16APR26_%j.err
#SBATCH --output=mmseqs2_80clustcov_16APR26_%j.out
#SBATCH --partition=k2-sandbox,k2-hipri
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --ntasks=16
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# ==============================================================================
# MMseqs2 all-vs-all search + representative clustering for network visualisation
# ==============================================================================
#
# This script:
#   1. Creates an MMseqs2 database from your protein FASTA
#   2. Clusters at 95% identity to get representative sequences
#   3. Runs an all-vs-all search of representatives against themselves
#   4. Outputs a TSV with pairwise similarities for network building
#   5. Outputs a cluster membership file so you can size nodes
#
# USAGE:
#   sbatch mmseqs2_defence_network.sh <input.fasta> <output_dir>
#
#   e.g. sbatch mmseqs2_defence_network.sh defence_proteins.fasta results/mmseqs2
# ==============================================================================

# --- Input validation ---
if [ $# -ne 2 ]; then
    echo "Usage: sbatch mmseqs2_80clustcov.sh <input.fasta> <output_dir>"
    exit 1
fi

INPUT_FASTA="$1"
OUTDIR="$2"

if [ ! -f "$INPUT_FASTA" ]; then
    echo "ERROR: Input FASTA not found: $INPUT_FASTA"
    exit 1
fi

# --- Setup ---
conda activate /mnt/scratch2/igfs-anaconda/conda-envs/mmseqs2_18.8/

mkdir -p "${OUTDIR}"
TMPDIR="${OUTDIR}/tmp"
mkdir -p "${TMPDIR}"

echo "============================================"
echo "MMseqs2 defence gene network pipeline"
echo "Input:    ${INPUT_FASTA}"
echo "Output:   ${OUTDIR}"
echo "Threads:  ${SLURM_NTASKS}"
echo "Started:  $(date)"
echo "============================================"

# --- Step 1: Create MMseqs2 database ---
echo "[$(date +%T)] Creating MMseqs2 database..."
mmseqs createdb "${INPUT_FASTA}" "${OUTDIR}/proteinDB"

# --- Step 2: Cluster at 95% identity + 80% coverage for representative sequences ---
echo "[$(date +%T)] Clustering at 95% identity to select representatives..."
mmseqs cluster \
    "${OUTDIR}/proteinDB" \
    "${OUTDIR}/clu95" \
    "${TMPDIR}" \
    --min-seq-id 0.95 \
    -c 0.8 \
    --cov-mode 3 \
    --threads "${SLURM_NTASKS}"

# Create TSV of cluster memberships (rep → member)
mmseqs createtsv \
    "${OUTDIR}/proteinDB" \
    "${OUTDIR}/proteinDB" \
    "${OUTDIR}/clu95" \
    "${OUTDIR}/cluster_membership.tsv"

# Extract representative sequences as a new database
mmseqs createsubdb \
    "${OUTDIR}/clu95" \
    "${OUTDIR}/proteinDB" \
    "${OUTDIR}/repDB"

# Also extract representative sequences as FASTA (useful for annotation)
mmseqs convert2fasta \
    "${OUTDIR}/repDB" \
    "${OUTDIR}/representative_seqs.fasta"

# Count representatives
N_REPS=$(grep -c "^>" "${OUTDIR}/representative_seqs.fasta")
N_TOTAL=$(grep -c "^>" "${INPUT_FASTA}")
echo "[$(date +%T)] Reduced ${N_TOTAL} proteins to ${N_REPS} representatives"

# --- Step 3: All-vs-all search of representatives ---
echo "[$(date +%T)] Running all-vs-all search of representatives..."
mmseqs search \
    "${OUTDIR}/repDB" \
    "${OUTDIR}/repDB" \
    "${OUTDIR}/result" \
    "${TMPDIR}" \
    -s 7.5 \
    --max-seqs 5000 \
    --threads "${SLURM_NTASKS}" \
    -e 1e-5

# --- Step 4: Convert results to TSV for network building ---
echo "[$(date +%T)] Converting results to TSV..."
mmseqs convertalis \
    "${OUTDIR}/repDB" \
    "${OUTDIR}/repDB" \
    "${OUTDIR}/result" \
    "${OUTDIR}/allvsall_results.tsv" \
    --format-output "query,target,pident,evalue,bits,qcov,tcov,alnlen"

# --- Step 5: Filter for network (80% coverage, 20% identity) ---
echo "[$(date +%T)] Filtering for network edges (pident>=20, coverage>=80)..."
awk -F'\t' '($3 >= 20) && ($6 >= 0.80) && ($7 >= 0.80) && ($1 != $2)' \
    "${OUTDIR}/allvsall_results.tsv" > "${OUTDIR}/network_edges.tsv"

# Add header to network edges file
sed -i '1i query\ttarget\tpident\tevalue\tbitscore\tqcov\ttcov\talnlen' \
    "${OUTDIR}/network_edges.tsv"

# --- Step 6: Create node size file (cluster membership counts) ---
echo "[$(date +%T)] Creating node size file from cluster memberships..."
awk -F'\t' '{print $1}' "${OUTDIR}/cluster_membership.tsv" \
    | sort | uniq -c | sort -rn \
    | awk '{print $2"\t"$1}' > "${OUTDIR}/node_sizes.tsv"
sed -i '1i representative\tcluster_size' "${OUTDIR}/node_sizes.tsv"

# --- Summary ---
N_EDGES=$(tail -n +2 "${OUTDIR}/network_edges.tsv" | wc -l)

echo ""
echo "============================================"
echo "Pipeline complete: $(date)"
echo "============================================"
echo "Total input proteins:    ${N_TOTAL}"
echo "Representative proteins: ${N_REPS}"
echo "Network edges:           ${N_EDGES}"
echo ""
echo "Output files:"
echo "  ${OUTDIR}/cluster_membership.tsv   - which proteins belong to which representative"
echo "  ${OUTDIR}/node_sizes.tsv           - representative ID + cluster size (for node sizing)"
echo "  ${OUTDIR}/representative_seqs.fasta - FASTA of representative sequences"
echo "  ${OUTDIR}/allvsall_results.tsv     - all pairwise hits (unfiltered)"
echo "  ${OUTDIR}/network_edges.tsv        - filtered edges ready for network visualisation"
echo ""
echo "Next steps:"
echo "  1. Load network_edges.tsv into Cytoscape/gephi/R (query-target = edge list)"
echo "  2. Join node_sizes.tsv to size nodes by cluster membership"
echo "  3. Colour nodes by defence system type from your annotations"
echo "============================================"

# --- Cleanup tmp ---
rm -rf "${TMPDIR}"
