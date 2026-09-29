#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Figure S4.3
# tblastn + clinker comparison of the strain prophages with Bonnie and Clyde.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/prophage_vs_phage_comparison.sh
# ======================================================================
# prophage_vs_phage_comparison.sh
#
# For each strain (C67, D32, D68):
#   1. Parse Bakta GFF to get coordinates of the prophage CDS features
#   2. Extract the prophage locus as a nucleotide FASTA (for blastn / clinker)
#   3. Extract prophage protein sequences from Bakta .faa (for tblastn)
#   4. Run tblastn: prophage proteins vs Bonnie and Clyde phage genomes
#
# Run from HPC:
#   bash prophage_vs_phage_comparison.sh
#
# Dependencies: samtools (faidx), seqkit, tblastn (BLAST+)
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
ANNOT_DIR="/mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations"
STRAIN_DIR="${ANNOT_DIR}/1_StrainGenomes"
PHAGE_DIR="${ANNOT_DIR}/2_PhageGenomes"
OUT_DIR="${ANNOT_DIR}/3_ProphageComparison"

mkdir -p "${OUT_DIR}"

# GFF annotations (Bakta output)
declare -A GFF
GFF["C67"]="${STRAIN_DIR}/C67_bakta/C67.gff"
GFF["D32"]="${STRAIN_DIR}/D32_bakta/D32.gff3"
GFF["D68"]="${STRAIN_DIR}/D68_bakta/D68.gff3"

# Genome nucleotide FASTAs (Bakta .fna)
declare -A FNA
FNA["C67"]="${STRAIN_DIR}/C67_bakta/C67.fna"
FNA["D32"]="${STRAIN_DIR}/D32_bakta/D32.fna"
FNA["D68"]="${STRAIN_DIR}/D68_bakta/D68.fna"

# Bakta protein FASTAs (.faa)
declare -A FAA
FAA["C67"]="${STRAIN_DIR}/C67_bakta/C67.faa"
FAA["D32"]="${STRAIN_DIR}/D32_bakta/D32.faa"
FAA["D68"]="${STRAIN_DIR}/D68_bakta/D68.faa"

# Phage genomes
BONNIE="${PHAGE_DIR}/BONNIE_prefixed.fasta"
CLYDE="${PHAGE_DIR}/CLYDE_prefixed.fasta"

# ── Prophage locus tags (from fig2_host_trajectories_FINALISED.R) ─────────────
declare -A PROPHAGE_IDS
PROPHAGE_IDS["C67"]="C67_prot_02009 C67_prot_00947 C67_prot_00948 C67_prot_00949 \
C67_prot_00950 C67_prot_00951 C67_prot_00952 C67_prot_00953 C67_prot_00954 \
C67_prot_00955 C67_prot_00956 C67_prot_00957 C67_prot_00958 C67_prot_00959 \
C67_prot_00960 C67_prot_00961 C67_prot_00962 C67_prot_00963 C67_prot_00964 \
C67_prot_00965"

PROPHAGE_IDS["D32"]="D32_prot_01483 D32_prot_01484 D32_prot_01485 D32_prot_01486 \
D32_prot_01487 D32_prot_01488 D32_prot_01489 D32_prot_01490 D32_prot_01491 \
D32_prot_01492 D32_prot_01493 D32_prot_01494 D32_prot_01495 D32_prot_01496 \
D32_prot_01497 D32_prot_01498 D32_prot_01499 D32_prot_01500 D32_prot_01501"

PROPHAGE_IDS["D68"]="D68_prot_01332 D68_prot_01333 D68_prot_01334 D68_prot_01335 \
D68_prot_01336 D68_prot_01337 D68_prot_01338 D68_prot_01339 D68_prot_01340 \
D68_prot_01341 D68_prot_01342 D68_prot_01343 D68_prot_01344 D68_prot_01345 \
D68_prot_01346 D68_prot_01347 D68_prot_01348 D68_prot_01349 D68_prot_01350"

# ── Step 1 & 2: Extract nucleotide and protein FASTAs per strain ──────────────
for STRAIN in C67 D32 D68; do
    echo ""
    echo "══ ${STRAIN} ═══════════════════════════════════════════════════════════"

    # Write locus tag list to file (used for GFF grep and seqkit)
    ID_FILE="${OUT_DIR}/${STRAIN}_prophage_ids.txt"
    echo "${PROPHAGE_IDS[$STRAIN]}" | tr ' ' '\n' > "${ID_FILE}"
    N_IDS=$(wc -l < "${ID_FILE}")
    echo "  Prophage genes: ${N_IDS}"

    # ── Extract matching CDS lines from GFF ──────────────────────────────────
    # Bakta uses locus_tag=<ID> in the attributes column (col 9)
    PROPH_GFF="${OUT_DIR}/${STRAIN}_prophage.gff3"
    grep $'\tCDS\t' "${GFF[$STRAIN]}" \
        | grep -Ff <(sed 's/^/locus_tag=/' "${ID_FILE}") \
        > "${PROPH_GFF}"

    N_CDS=$(wc -l < "${PROPH_GFF}")
    echo "  CDS features found in GFF: ${N_CDS}"

    if [[ "${N_CDS}" -eq 0 ]]; then
        echo "  ERROR: no CDS lines matched — check locus tag format in GFF"
        continue
    fi

    # ── Get locus span: contig, min start, max end (GFF cols 1, 4, 5) ────────
    # GFF3 coordinates are 1-based, samtools faidx expects 1-based region syntax
    read -r CONTIG LOCUS_START LOCUS_END < <(
        awk 'BEGIN{min=999999999; max=0; seq=""}
             {if(NR==1) seq=$1;
              if($4+0 < min) min=$4+0;
              if($5+0 > max) max=$5+0}
             END{print seq, min, max}' "${PROPH_GFF}"
    )
    LOCUS_LEN=$(( LOCUS_END - LOCUS_START + 1 ))
    echo "  Prophage locus: ${CONTIG}:${LOCUS_START}-${LOCUS_END} (${LOCUS_LEN} bp)"

    # ── Extract nucleotide FASTA ──────────────────────────────────────────────
    OUT_FNA="${OUT_DIR}/${STRAIN}_prophage_locus.fna"
    samtools faidx "${FNA[$STRAIN]}" "${CONTIG}:${LOCUS_START}-${LOCUS_END}" \
        | sed "1s/.*/>${STRAIN}_prophage  ${CONTIG}:${LOCUS_START}-${LOCUS_END}  ${LOCUS_LEN}bp/" \
        > "${OUT_FNA}"
    echo "  Nucleotide FASTA → ${OUT_FNA}"

    # ── Extract protein FASTA ─────────────────────────────────────────────────
    OUT_FAA="${OUT_DIR}/${STRAIN}_prophage_proteins.faa"
    samtools faidx "${FAA[$STRAIN]}"
    xargs samtools faidx "${FAA[$STRAIN]}" < "${ID_FILE}" > "${OUT_FAA}"
    N_PROTS=$(grep -c "^>" "${OUT_FAA}")
    echo "  Protein FASTA → ${OUT_FAA} (${N_PROTS} sequences)"

    if [[ "${N_PROTS}" -lt "${N_IDS}" ]]; then
        echo "  WARNING: only ${N_PROTS}/${N_IDS} proteins found in FAA — check IDs"
    fi
done

# ── Step 3: tblastn — prophage proteins vs phage genomes ─────────────────────
echo ""
echo "══ tblastn: prophage proteins vs Bonnie and Clyde ══════════════════════"

# Build BLAST databases for phage genomes (only once)
for PHAGE_GENOME in "${BONNIE}" "${CLYDE}"; do
    if [[ ! -f "${PHAGE_GENOME}.nhr" ]]; then
        echo "  Making BLAST db: $(basename ${PHAGE_GENOME})"
        makeblastdb -in "${PHAGE_GENOME}" -dbtype nucl -parse_seqids
    fi
done

TBLASTN_HEADER="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"

for STRAIN in C67 D32 D68; do
    QUERY="${OUT_DIR}/${STRAIN}_prophage_proteins.faa"
    [[ -f "${QUERY}" ]] || { echo "  SKIP ${STRAIN}: protein FASTA missing"; continue; }

    for PHAGE_NAME in BONNIE CLYDE; do
        case "${PHAGE_NAME}" in
            BONNIE) SUBJECT="${BONNIE}" ;;
            CLYDE)  SUBJECT="${CLYDE}"  ;;
        esac

        OUT_BLAST="${OUT_DIR}/${STRAIN}_prophage_vs_${PHAGE_NAME}.tsv"
        echo "  ${STRAIN} prophage → ${PHAGE_NAME}..."

        # Write header then results
        echo -e "${TBLASTN_HEADER}" > "${OUT_BLAST}"
        tblastn \
            -query   "${QUERY}" \
            -db      "${SUBJECT}" \
            -evalue  1e-5 \
            -num_threads 4 \
            -outfmt  "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            >> "${OUT_BLAST}"

        N_HITS=$(( $(wc -l < "${OUT_BLAST}") - 1 ))
        echo "    ${N_HITS} hits → ${OUT_BLAST}"
    done
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "══ Output files ════════════════════════════════════════════════════════"
echo "  Nucleotide loci (for clinker / blastn):"
ls "${OUT_DIR}"/*_prophage_locus.fna 2>/dev/null | sed 's/^/    /'
echo ""
echo "  Protein FASTAs (for tblastn):"
ls "${OUT_DIR}"/*_prophage_proteins.faa 2>/dev/null | sed 's/^/    /'
echo ""
echo "  tblastn results:"
ls "${OUT_DIR}"/*_prophage_vs_*.tsv 2>/dev/null | sed 's/^/    /'
echo ""
echo "  For clinker: pass the three *_prophage_locus.fna files + BONNIE/CLYDE fastas:"
echo "    clinker ${OUT_DIR}/*_prophage_locus.fna ${BONNIE} ${CLYDE} -p clinker_output.html"
echo ""
echo "Done."

