#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Figures S4.4-S4.6
# Pulls prophage loci (+5 kb flanks) out of the per-sample coverage bedgraphs.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/4_SuccessfulRuns/Coverages/extract_prophage_bedgraph.sh
# ======================================================================
# extract_prophage_bedgraph.sh
# Extracts prophage locus + 5 kb flanking from all bedgraph files.
# Outputs one tiny TSV per sample + a normalization factors table.

set -euo pipefail

COV_DIR="/mnt/scratch2/users/40204129/3_Transcriptomics/4_SuccessfulRuns/Coverages"
OUT_DIR="/mnt/scratch2/users/40204129/3_Transcriptomics/4_SuccessfulRuns/prophage_coverage"
mkdir -p "${OUT_DIR}"

# Prophage coords (exact) + 5 kb flanking either side
# Format: CONTIG  FLANK_START  PROPH_START  PROPH_END  FLANK_END
declare -A CONTIG FS PS PE FE
CONTIG[C67]="contig_4";  FS[C67]=50949; PS[C67]=55949; PE[C67]=69028; FE[C67]=74028
CONTIG[D32]="contig_10"; FS[D32]=291;   PS[D32]=5291;  PE[D32]=18370; FE[D32]=23370
CONTIG[D68]="contig_11"; FS[D68]=291;   PS[D68]=5291;  PE[D68]=18370; FE[D68]=23370

# Header for normalization file (total mapped bases per sample, for RPM scaling)
NORM="${OUT_DIR}/normalization_factors.tsv"
echo -e "sample\ttimepoint\tstrain\tcondition\treplicate\ttotal_bases_mapped" > "${NORM}"

for BG in "${COV_DIR}"/*.bedgraph; do
    FNAME=$(basename "${BG}" __coverage.bedgraph)

    # Parse filename: {timepoint}-{strain}-{condition}-{replicate}
    IFS='-' read -r TP STRAIN COND REP <<< "${FNAME}"

    if [[ -z "${CONTIG[$STRAIN]+x}" ]]; then
        echo "  SKIP ${FNAME}: unrecognised strain '${STRAIN}'"
        continue
    fi

    C="${CONTIG[$STRAIN]}"
    RSTART="${FS[$STRAIN]}"   # flanked region start
    REND="${FE[$STRAIN]}"     # flanked region end

    OUT="${OUT_DIR}/${FNAME}.tsv"

    # Single awk pass:
    #   - accumulates total_bases_mapped across the whole file (for normalisation)
    #   - prints region-overlapping intervals to the output TSV
    awk -v contig="$C" -v rstart="$RSTART" -v rend="$REND" \
        -v tp="$TP" -v strain="$STRAIN" -v cond="$COND" -v rep="$REP" \
        -v outfile="${OUT}" -v sample="$FNAME" -v normfile="${NORM}" \
    'BEGIN { total = 0; OFS = "\t" }
     { total += $4 * ($3 - $2) }
     $1 == contig && $3 > rstart && $2 < rend {
         # clip interval to the requested window
         s = ($2 < rstart) ? rstart : $2
         e = ($3 > rend)   ? rend   : $3
         print tp, strain, cond, rep, s, e, $4 >> outfile
     }
     END {
         print sample, tp, strain, cond, rep, total >> normfile
     }' "${BG}"

    echo "  ${FNAME}: $(wc -l < "${OUT}") intervals extracted"
done

echo ""
echo "Done. Download the whole folder:"
echo "  scp -r ${OUT_DIR} your_local_path/"
